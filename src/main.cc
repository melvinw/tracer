#include <cassert>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>
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
#include <glm/gtx/rotate_vector.hpp>
#include <json.hpp>

#include "gltf_loader.h"
#include "gpl/solar_position.h"
#include "types.h"

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
  Scalar absorbed_flux;   // Watts / m^2
  Scalar scattered_flux;  // Watts / m^2
};

std::optional<std::pair<struct tm, int>> parse_time(
    const std::string &time_string) {
  int tz;
  struct tm tv;
  char *end = strptime(time_string.c_str(), "%Y-%m-%dT%H:%M:%S", &tv);
  struct tm tv2;
  if (end == nullptr) {
    std::cerr << "Error in parsing time " << time_string << std::endl;
    return std::nullopt;
  }
  if (*end == '-' || *end == '+') {
    end = strptime(end + 1, "%H:%M", &tv2);
    assert(end != nullptr);
    tz = (*end == '-') ? -tv2.tm_hour : tv2.tm_hour;
  }
  return std::make_pair(tv, tz);
}

int main(int argc, char **argv) {
  // TODO: use a real logging library
  CLI::App app{"Tracer"};
  bool debug = getenv("DEBUG") != nullptr;

  std::string gltf_path;
  app.add_option("--gltf_path", gltf_path, "Path to GLTF")->required();

  double latitude = 0.0f;
  double longitude = 0.0f;
  app.add_option("--lat", latitude, "Latitude")->required();
  app.add_option("--long", longitude, "Longitude")->required();

  std::string time_string;
  CLI::Option *time_option = app.add_option(
      "--time", time_string,
      "Time of day - should be in %Y-%m-%dT%H:%M%:S[-/+%H:%M] format");

  CLI11_PARSE(app, argc, argv);

  struct tm tv;
  tv.tm_hour = 12;
  int tz = 0;
  if (time_option->count() > 0) {
    auto time = parse_time(time_string);
    assert(time.has_value());
    std::tie(tv, tz) = time.value();
  }

  if (debug) {
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
  float zenith_angle = sunAngle(latitude, longitude, tz, tv);
  glm::vec3 sun_center_glm = glm::rotateY(
      glm::vec3(sun_center[0], sun_center[1], sun_center[2]), zenith_angle);
  sun_center = Vec3(sun_center_glm[0], sun_center_glm[1], sun_center_glm[2]);
  Vec3 sun_norm = bvh::normalize(sun_center - bsphere.origin);
  Scalar d = -bvh::dot(sun_center, sun_norm);
  Scalar sun_flux = 500 * std::cos(zenith_angle);  // W/m^2
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
