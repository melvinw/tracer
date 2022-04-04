#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include<math.h>
#include<time.h>
#include <bvh/bounding_box.hpp>
#include <bvh/bvh.hpp>
#include <bvh/heuristic_primitive_splitter.hpp>
#include <bvh/leaf_collapser.hpp>
#include <bvh/node_layout_optimizer.hpp>
#include <bvh/parallel_reinsertion_optimizer.hpp>
#include <bvh/primitive_intersectors.hpp>
#include <bvh/single_ray_traverser.hpp>
#include <bvh/sphere.hpp>
#include <bvh/sweep_sah_builder.hpp>
#include <bvh/triangle.hpp>
#include <bvh/vector.hpp>
#include <json.hpp>

#include "gltf_loader.h"
#include "types.h"

#include <glm/gtx/rotate_vector.hpp>
using Scalar = double;
using Vec3 = bvh::Vector3<Scalar>;
using Triangle = bvh::Triangle<Scalar>;
using Bvh = bvh::Bvh<Scalar>;
using BoundingBox = bvh::BoundingBox<Scalar>;
using Sphere = bvh::Sphere<Scalar>;

struct Stats {
  uint32_t hits;
  uint32_t self_hits;
  uint32_t misses;
  Scalar absorbed_flux;  // Watts / m^2
  Scalar scattered_flux; // Watts / m^2
};

//Taken from Helios
unsigned int julianDay(struct tm time_struct){
  unsigned int day = time_struct.tm_mday;
  unsigned int month = time_struct.tm_mon + 1; //Code below assumes months from 1-12 
  unsigned int year = time_struct.tm_year+1900; //Time struct contains years relative to the year 1900
  int skips_leap[] = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};
  int skips_nonleap[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  int* skips;
  
  if( (year-2000)%4 == 0 ){  //leap year
    skips=skips_leap;
  }else{                 //non-leap year
    skips=skips_nonleap;
  }
  
  return skips[month-1]+day;
}

double acos_safe( double x ){
  if (x < -1.0) x = -1.0 ;
  else if (x > 1.0) x = 1.0 ;
  return acosf(x) ;
}

double asin_safe( double x ){
  if (x < -1.0) x = -1.0 ;
  else if (x > 1.0) x = 1.0 ;
  return asinf(x) ;
}
//Returns the Zenith angle of the sun. Assumes phi = 0 for now
//From https://github.com/PlantSimulationLab/Helios/blob/master/plugins/solarposition/src/SolarPosition.cpp
double sunAngle(double latitude, double longitude, int UTC, struct tm time_struct){
  unsigned int jd = julianDay(time_struct);

  double rad=M_PI/180.f;
  double Gamma = 2.f*M_PI*(float(jd-1))/365.f;  
  //solar declination angle (Iqbal Eq. 1.3.1 after Spencer)
  double delta = 0.006918f - 0.399912f*cos(Gamma) + 0.070257f*sin(Gamma) 
    - 0.006758f*cos(2.f*Gamma) + 0.000907f*sin(2.f*Gamma) 
    - 0.002697f*cos(3.f*Gamma) + 0.00148f*sin(3.f*Gamma);
  double EoT = 229.18f*(0.000075f + 0.001868f*cos(Gamma) - 0.032077f*sin(Gamma)
		 - 0.014615f*cos(2.f*Gamma) - 0.04089f*sin(2.f*Gamma));
  double time_dec=time_struct.tm_hour+time_struct.tm_min/60.f;  //(hours) 

  int LSTM=15.f*float(UTC); //degrees

  double TC=4.f*(LSTM-longitude)+EoT; //minutes
  double LST=time_dec+TC/60.f; //hours

  double h=(LST-12.f)*15.f*rad; //hour angle (rad)
    
  //solar zenith angle
  double theta = asin_safe( sin(latitude*rad)*sin(delta) + cos(latitude*rad)*cos(delta)*cos(h) ); //(rad)

  assert( theta>-0.5f*M_PI && theta<0.5f*M_PI );

  //solar elevation angle
  double phi = acos_safe( (sin(delta) - sin(theta)*sin(latitude*rad))/(cos(theta)*cos(latitude*rad)));

  if( LST>12.f ){
    phi=2.f*M_PI-phi;
  }

  assert( phi>0 && phi<2.f*M_PI );

  return (90 * rad) - theta; //Ignoring phi for now. 90 - theta to return zenith angle
  
}

void sunAngleTest(){
  //Test taken from Helios
  double test_latitude = 40.1250;
  double test_longitude = 105.2369;
  struct tm* test_tm;
  time_t test_time;
  time(&test_time);
  test_tm = localtime(&test_time);
  test_tm->tm_mday = 1;
  test_tm->tm_mon = 0;
  test_tm->tm_year = 100;
  test_tm->tm_hour = 10;
  test_tm->tm_min = 30;
  test_tm->tm_sec = 0;
  int UTC = 7;
  std::cout << "Azimuthal angle for 10AM in CO time zone is: " << sunAngle(test_latitude, test_longitude, UTC, *test_tm) << std::endl;
}
int main(int argc, char **argv) {
  // TODO: use a real logging library
  bool debug = getenv("DEBUG") != nullptr;
  

  sunAngleTest();
  
  
  assert(argc >= 2);
  std::string gltf_path(argv[1]);
  double latitude = std::stof(argv[2]);
  double longitude = std::stof(argv[3]);
  int UTC = std::stoi(argv[4]);
  struct tm time_value;
  char* err = strptime(argv[5], "%Y-%m-%d %H:%M:%S", &time_value);
  if(err == NULL || *err != '\0'){
    std::cerr << "Error in parsing time " << argv[5] << std::endl;
  }
  if (debug) {
    const int buf_sz = 1024;
    char buf[buf_sz];
    strftime(buf, buf_sz, "%c", &time_value);
    std::cerr << "Time is: " << buf << std::endl;
    std::cerr << "loading " << gltf_path << std::endl;
  }
  
  auto [meshes, mesh_instances] = LoadMeshesFromGLTF(gltf_path);

  std::vector<Triangle> triangles;
  std::map<int, std::pair<size_t, size_t>> instance_triangles;
  for (size_t idx = 0; idx < mesh_instances.size(); idx++) {
    const auto &instance = mesh_instances[idx];
    const auto &mesh = meshes[instance.mesh_index];
    size_t start = triangles.size();
    for (size_t i = 0; i < mesh.indices.size(); i += 3) {
      glm::vec4 v0 = mesh.vertices[mesh.indices[i]].position;
      glm::vec4 v1 = mesh.vertices[mesh.indices[i + 1]].position;
      glm::vec4 v2 = mesh.vertices[mesh.indices[i + 2]].position;
      if (instance.model_matrix.has_value()) {
        v0 = instance.model_matrix.value() * v0;
        v1 = instance.model_matrix.value() * v1;
        v2 = instance.model_matrix.value() * v2;
      }
      triangles.push_back(Triangle(Vec3(v0.x, v0.y, v0.z),
                                   Vec3(v1.x, v1.y, v1.z),
                                   Vec3(v2.x, v2.y, v2.z)));
    }
    size_t end = triangles.size() - 1;
    instance_triangles[idx] = {start, end};
  }

  if (debug) {
    std::cerr << "found " << triangles.size() << " triangles comprising "
              << mesh_instances.size() << " models. constructing BVH"
              << std::endl;
  }
  Bvh bvh;
  auto [bboxes, centers] = bvh::compute_bounding_boxes_and_centers(
      triangles.data(), triangles.size());
  auto global_bbox =
      bvh::compute_bounding_boxes_union(bboxes.get(), triangles.size());

  bvh::SweepSahBuilder<Bvh> builder(bvh);
  builder.build(global_bbox, bboxes.get(), centers.get(), triangles.size());

  bvh::ParallelReinsertionOptimizer<Bvh> reinsertion_optimizer(bvh);
  reinsertion_optimizer.optimize();

  bvh::NodeLayoutOptimizer layout_optimizer(bvh);
  layout_optimizer.optimize();

  Scalar radius = bvh::length(global_bbox.diagonal()) / 2.0;
  Sphere bsphere(global_bbox.center(), radius);
  if (debug) {
    std::cerr << "bounding sphere at (" << bsphere.origin[0] << ", "
              << bsphere.origin[1] << ", " << bsphere.origin[2]
              << ") w/ radius " << bsphere.radius << std::endl;
  }

  // TODO: Sun is fixed at high noon right now. Make this configurable.
  Vec3 sun_center = bsphere.origin + Vec3(0, 0, bsphere.radius);
  float zenith_angle = sunAngle(latitude, longitude, UTC, time_value);
  glm::vec3 sun_center_glm = glm::rotateY(glm::vec3(sun_center[0], sun_center[1], sun_center[2]), zenith_angle);
  sun_center = Vec3(sun_center_glm[0], sun_center_glm[1], sun_center_glm[2]);
  Vec3 sun_norm = bvh::normalize(sun_center - bsphere.origin);
  Scalar d = -bvh::dot(sun_center, sun_norm);
  Scalar sun_flux = 500 * std::cos(zenith_angle); // W/m^2
  constexpr Scalar absorb_factor = Scalar(3) / Scalar(4);
  // TODO scatter factor should come from material properties of the underlying
  // mesh
  constexpr Scalar scatter_factor = Scalar(1) / Scalar(10);
  constexpr size_t rays_per_triangle = 1;
  if (debug) {
    std::cerr << "Zenith angle is: " << zenith_angle << std::endl;
    std::cerr << "sun disk at (" << sun_center[0] << ", " << sun_center[1]
              << ", " << sun_center[2] << ")" << std::endl;
  }
  // First pass to accumulate light on each mesh
  bvh::ClosestPrimitiveIntersector<Bvh, Triangle> primitive_intersector(
      bvh, triangles.data());
  bvh::SingleRayTraverser<Bvh> traverser(bvh);
  std::vector<Stats> stats(mesh_instances.size());
  for (const auto &[idx, interval] : instance_triangles) {
    const auto &instance = mesh_instances[idx];
    const auto [start, end] = interval;
    if (debug) {
      std::cerr << "launching " << end - start + 1 << " rays from "
                << instance.name << "(" << start << ", " << end << ")"
                << std::endl;
    }

    Stats &s = stats[idx];
    for (size_t i = start; i <= end; i++) {
      const auto &t = triangles[i];

      for (size_t j = 0; j < rays_per_triangle; j++) {
        // offset ray origin by scaled down normal to avoid self intersections.
        // TODO sample ray origins from surface instead of casting from center
        auto origin = t.center() + t.n * .000000001;
        auto dir = origin - sun_norm * (bvh::dot(sun_norm, origin) + d);
        bvh::Ray<Scalar> ray(origin, bvh::normalize(dir), 0.000001,
                             2.0 * bsphere.radius);

        auto hit = traverser.traverse(ray, primitive_intersector);
        if (hit.has_value()) {
          auto tidx = hit->primitive_index;
          if (tidx >= start && tidx <= end) {
            s.self_hits++;
          } else {
            s.hits++;
          }
        } else {
          s.misses++;
          // for details on this math see pp14 of
          // https://www.sciencedirect.com/science/article/pii/S0304380017304842
          Scalar absorbed = absorb_factor * sun_flux *
                            glm::abs(bvh::dot(t.n, sun_norm)) /
                            Scalar(rays_per_triangle);
          s.absorbed_flux += absorbed;
          s.scattered_flux += scatter_factor * absorbed;
        }
      }
    }
  }

  // TODO: Second pass to scatter some portion of light form each surface

  nlohmann::json output;
  for (size_t i = 0; i < mesh_instances.size(); i++) {
    const auto &instance = mesh_instances[i];
    const auto &s = stats[i];
    output[instance.name] = {
        {"obstructed_rays", s.hits + s.self_hits},
        {"unobstructed_rays", s.misses},
        {"absorbed_flux", s.absorbed_flux},
        {"scattered_flux", s.scattered_flux},
    };
  }
  std::cout << output.dump(4) << std::endl;

  return 0;
}
