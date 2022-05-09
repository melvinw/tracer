#include <atomic>
#include <cassert>
#include <ctime>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <thread>
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
#include <thread_pool.hpp>

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
  Stats()
      : hits(0),
        misses(0),
        absorbed_flux(Scalar(0.0f)),
        scattered_flux(Scalar(0.0f)) {}

  Stats(const Stats &other)
      : hits(other.hits.load()),
        misses(other.misses.load()),
        absorbed_flux(other.absorbed_flux.load()),
        scattered_flux(other.absorbed_flux.load()) {}

  Stats(Stats &&other)
      : hits(other.hits.load()),
        misses(other.misses.load()),
        absorbed_flux(other.absorbed_flux.load()),
        scattered_flux(other.absorbed_flux.load()) {}
  Stats &operator=(Stats &&other) {
    hits = other.hits.load();
    misses = other.misses.load();
    absorbed_flux = other.absorbed_flux.load();
    scattered_flux = other.absorbed_flux.load();
    return *this;
  }
  void Merge(const Stats &s) {
    hits = hits.load() + s.hits.load();
    misses = misses.load() + s.misses.load();
    absorbed_flux = absorbed_flux.load() + s.absorbed_flux.load();
    scattered_flux = scattered_flux.load() + s.scattered_flux.load();
  }

  std::atomic<uint32_t> hits;
  std::atomic<uint32_t> misses;
  std::atomic<Scalar> absorbed_flux;   // Watts / m^2
  std::atomic<Scalar> scattered_flux;  // Watts / m^2
};

struct AnnotatedTriangle {
  struct Intersection {
    Scalar t;

    // Required member: returns the distance along the ray
    Scalar distance() const { return t; }
  };

  using ScalarType = Scalar;
  using IntersectionType = Intersection;

  AnnotatedTriangle() = default;
  AnnotatedTriangle(size_t idx, const Triangle &t)
      : mesh_idx(idx), triangle(t), stats(){};

  Vec3 center() const { return triangle.center(); }

  BoundingBox bounding_box() const { return triangle.bounding_box(); }

  std::optional<Intersection> intersect(const bvh::Ray<Scalar> &ray) const {
    if (auto ret = triangle.intersect(ray); ret.has_value()) {
      return Intersection{.t = ret.value().t};
    }
    return std::nullopt;
  }

  size_t mesh_idx;
  Triangle triangle;
  Stats stats;
};

// TODO: move this to types.h
struct Plane {
  Plane() : a(), b(), c(), d() {}

  // normal is assumed to be unit lenght.
  Plane(const Vec3 &point, const Vec3 &normal)
      : a(normal[0]), b(normal[1]), c(normal[2]), d(-bvh::dot(point, normal)) {}

  Vec3 ProjectPoint(const Vec3 &point) {
    Vec3 normal(a, b, c);
    Scalar signed_dist = bvh::dot(normal, point) + d;
    return point - normal * signed_dist;
  }

  double a, b, c, d;
};

static std::pair<struct tm, int> parse_time(const std::string &time_string) {
  int tz;
  struct tm tv;
  char *end = strptime(time_string.c_str(), "%Y-%m-%dT%H:%M:%S", &tv);
  struct tm tv2;
  assert(end != nullptr);
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

  std::string gltf_path;
  app.add_option("--gltf_path", gltf_path, "Path to GLTF")->required();

  double latitude = 0.0f;
  double longitude = 0.0f;
  app.add_option("--lat", latitude, "Latitude")->required();
  app.add_option("--long", longitude, "Longitude")->required();

  size_t rays_per_triangle = 0;
  app.add_option("--rays-per-triangle", rays_per_triangle,
                 "Number of rays to cast from each triangle in the mesh")
      ->default_val(1);
  assert(rays_per_triangle != 0);

  unsigned int num_workers = 1;
  app.add_option("--num-workers", num_workers,
                 "Number of workers to use when casting rays.")
      ->default_val(std::thread::hardware_concurrency());
  assert(num_workers != 0);

  std::string start_time_string;
  CLI::Option *start_time_option = app.add_option(
      "--time", start_time_string,
      "Time of day - should be in %Y-%m-%dT%H:%M:%S[-/+%H:%M] format");

  std::string end_time_string;
  CLI::Option *end_time_option = app.add_option(
      "--end-time", end_time_string,
      "End time - should be in %Y-%m-%dT%H:%M:%S[-/+%H:%M] format");

  std::string increment_time_string;
  CLI::Option *increment_time_option =
      app.add_option("--increment-time", increment_time_string,
                     "Increment time - should be in %Y-%m-%dT%H:%M:%Sformat");
  bool debug = false;
  app.add_flag("--verbose", debug, "Be verbose.");

  CLI11_PARSE(app, argc, argv);

  struct tm start_tv = {};
  start_tv.tm_mon = 1;
  start_tv.tm_mday = 1;
  start_tv.tm_hour = 12;
  int start_tz = 0;
  if (start_time_option->count() > 0) {
    std::tie(start_tv, start_tz) = parse_time(start_time_string);
  }

  struct tm end_tv = start_tv;
  int end_tz = end_tz;
  if (end_time_option->count() > 0) {
    std::tie(end_tv, end_tz) = parse_time(end_time_string);
  }

  struct tm increment_tv;
  increment_tv.tm_mon = 0;
  increment_tv.tm_mday = 0;
  increment_tv.tm_hour = 1;  // Increments by hour by default
  increment_tv.tm_min = 0;
  increment_tv.tm_sec = 0;
  int increment_tz;
  if (increment_time_option->count() > 0) {
    std::tie(increment_tv, increment_tz) =
        parse_time(increment_time_string);  // Increment_tz is irrelevant
  }

  if (debug) {
    std::cerr << "loading " << gltf_path << std::endl;
  }

  auto [meshes, mesh_instances] = LoadMeshesFromGLTF(gltf_path);

  std::vector<AnnotatedTriangle> triangles;
  for (size_t idx = 0; idx < mesh_instances.size(); idx++) {
    const auto &instance = mesh_instances[idx];
    const auto &mesh = meshes[instance.mesh_index];
    for (size_t i = 0; i < mesh.indices.size(); i += 3) {
      glm::vec4 v0 = mesh.vertices[mesh.indices[i]].position;
      glm::vec4 v1 = mesh.vertices[mesh.indices[i + 1]].position;
      glm::vec4 v2 = mesh.vertices[mesh.indices[i + 2]].position;
      if (instance.model_matrix.has_value()) {
        v0 = instance.model_matrix.value() * v0;
        v1 = instance.model_matrix.value() * v1;
        v2 = instance.model_matrix.value() * v2;
      }
      AnnotatedTriangle t(
          idx, Triangle(Vec3(v0.x, v0.y, v0.z), Vec3(v1.x, v1.y, v1.z),
                        Vec3(v2.x, v2.y, v2.z)));
      triangles.emplace_back(std::move(t));
    }
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

  Vec3 sun_center = bsphere.origin + Vec3(0, 0, bsphere.radius);
  struct tm current_tv = start_tv;

  // mktime takes care of value overflows for struct tm values and modifies the
  // struct tm being passed in
  while (mktime(&current_tv) <= mktime(&end_tv)) {
    float zenith_angle = sunAngle(latitude, longitude, start_tz, current_tv);
    glm::vec3 sun_center_glm = glm::rotateY(
        glm::vec3(sun_center[0], sun_center[1], sun_center[2]), zenith_angle);
    sun_center = Vec3(sun_center_glm[0], sun_center_glm[1], sun_center_glm[2]);
    Vec3 sun_norm = bvh::normalize(sun_center - bsphere.origin);
    Plane sun_plane(sun_center, sun_norm);
    Scalar sun_flux = 1362 / std::cos(zenith_angle);  // W/m^2

    // TODO scatter factor should come from material properties of the
    // underlying mesh
    constexpr Scalar scatter_factor = Scalar(1) / Scalar(10);
    constexpr Scalar absorb_factor = Scalar(3) / Scalar(4);

    if (debug) {
      std::cerr << "Zenith angle is: " << zenith_angle << std::endl;
      std::cerr << "sun disk at (" << sun_center[0] << ", " << sun_center[1]
                << ", " << sun_center[2] << ")" << std::endl;
    }

    thread_pool tpool(num_workers);

    // First pass to accumulate light on each mesh
    bvh::AnyPrimitiveIntersector<Bvh, AnnotatedTriangle> intersector(
        bvh, triangles.data());
    bvh::SingleRayTraverser<Bvh> traverser(bvh);
    tpool.parallelize_loop(
        0ull, triangles.size(), [&](const size_t &lo, const size_t &hi) {
          for (size_t i = lo; i < hi; i++) {
            auto &t = triangles[i];

            for (size_t j = 0; j < rays_per_triangle; j++) {
              // offset ray origin by scaled down normal to avoid self
              // intersections.
              // TODO sample ray origins from surface instead of casting from
              // center
              auto origin = t.triangle.center() + t.triangle.n * .000000001;
              auto dir = origin - sun_plane.ProjectPoint(origin);
              bvh::Ray<Scalar> ray(origin, bvh::normalize(dir), 0.000001,
                                   2.0 * bsphere.radius);

              auto hit = traverser.traverse(ray, intersector);
              if (hit.has_value()) {
                t.stats.hits++;
              } else {
                t.stats.misses++;
                // for details on this math see pp14 of
                // https://www.sciencedirect.com/science/article/pii/S0304380017304842
                Scalar absorbed = absorb_factor * sun_flux *
                                  glm::abs(bvh::dot(t.triangle.n, sun_norm)) /
                                  Scalar(rays_per_triangle);
                // C++17 doesn't have the std::atomic specializations for
                // floating point types yet, so we just repeatedly CAS here
                // until we succeed.
                for (Scalar af = t.stats.absorbed_flux;
                     !t.stats.absorbed_flux.compare_exchange_strong(
                         af, af + absorbed);) {
                  ;
                }
                for (Scalar sf = t.stats.scattered_flux;
                     !t.stats.scattered_flux.compare_exchange_strong(
                         sf, sf + scatter_factor * absorbed);) {
                  ;
                }
              }
            }
          }
        });
    nlohmann::json output;
    std::vector<Stats> stats(mesh_instances.size());
    for (const auto &t : triangles) {
      stats[t.mesh_idx].Merge(t.stats);
    }
    for (size_t i = 0; i < mesh_instances.size(); i++) {
      const auto &instance = mesh_instances[i];
      const auto &s = stats[i];
      output[instance.name] = {
          {"obstructed_rays", s.hits.load()},
          {"unobstructed_rays", s.misses.load()},
          {"absorbed_flux", s.absorbed_flux.load()},
          {"scattered_flux", s.scattered_flux.load()},
      };
    }
    for (auto &t : triangles) {
      t.stats = Stats();
    }
    std::cout << output.dump(4) << std::endl;

    // TODO: Second pass to scatter some portion of light form each surface
    current_tv.tm_mon += increment_tv.tm_mon;
    current_tv.tm_mday += increment_tv.tm_mday;
    current_tv.tm_hour += increment_tv.tm_hour;
    current_tv.tm_min += increment_tv.tm_min;
    current_tv.tm_sec += increment_tv.tm_sec;
  }

  return 0;
}
