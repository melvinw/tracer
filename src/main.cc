#include <cassert>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

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

#include "gltf_loader.h"
#include "types.h"

using Scalar = double;
using Vec3 = bvh::Vector3<Scalar>;
using Triangle = bvh::Triangle<Scalar>;
using Bvh = bvh::Bvh<Scalar>;
using BoundingBox = bvh::BoundingBox<Scalar>;
using Sphere = bvh::Sphere<Scalar>;

int main(int argc, char **argv) {
  assert(argc == 2);
  std::string gltf_path(argv[1]);
  std::cerr << "loading " << gltf_path << std::endl;
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

  std::cerr << "found " << mesh_instances.size() << " models consisting of "
            << triangles.size() << " triangles. constructing BVH" << std::endl;
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
  std::cerr << "bounding sphere at (" << bsphere.origin[0] << ", "
            << bsphere.origin[1] << ", " << bsphere.origin[2] << ") w/ radius "
            << bsphere.radius << std::endl;

  // TODO: Sun is fixed at high noon right now. Make this configurable.
  Vec3 sun_center = bsphere.origin + Vec3(0, 0, bsphere.radius);
  std::cerr << "sun disk at (" << sun_center[0] << ", " << sun_center[1] << ", "
            << sun_center[2] << ")" << std::endl;

  // First pass to accumulate light on each mesh
  for (const auto &[idx, interval] : instance_triangles) {
    const auto &instance = mesh_instances[idx];
    const auto [start, end] = interval;
    std::cerr << "launching " << end - start + 1 << " rays from "
              << instance.name << "(" << start << ", " << end << ")"
              << std::endl;

    size_t hits = 0;
    size_t self_hits = 0;
    size_t misses = 0;
    for (size_t i = start; i <= end; i++) {
      const auto &t = triangles[i];

      // offset ray origin by scaled down normal to avoid self intersections.
      auto origin = t.center() + t.n * .000000001;
      // TODO: cast to closest point on disk instead of center
      bvh::Ray<Scalar> ray(origin, bvh::normalize(sun_center - origin),
                           0.000001, 2.0 * bsphere.radius);
      bvh::ClosestPrimitiveIntersector<Bvh, Triangle> primitive_intersector(
          bvh, triangles.data());
      bvh::SingleRayTraverser<Bvh> traverser(bvh);

      if (auto hit = traverser.traverse(ray, primitive_intersector);
          hit.has_value()) {
        auto tidx = hit->primitive_index;
        if (tidx >= start && tidx <= end) {
          self_hits++;
        } else {
          hits++;
        }
      } else {
        misses++;
      }
    }
    std::cerr << "\thits: " << hits << ", self_hits: " << self_hits
              << ", misses: " << misses << std::endl;
  }

  // TODO: Second pass to scatter some portion of light form each surface

  // TODO: output copy of gltf from user w/ radiance annotations

  return 0;
}
