#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <bvh/bounding_box.hpp>
#include <bvh/bvh.hpp>
#include <bvh/primitive_intersectors.hpp>
#include <bvh/single_ray_traverser.hpp>
#include <bvh/sphere.hpp>
#include <bvh/sweep_sah_builder.hpp>
#include <bvh/triangle.hpp>
#include <bvh/vector.hpp>
#include <tiny_gltf.h>

using Scalar = float;
using Vec3 = bvh::Vector3<Scalar>;
using Triangle = bvh::Triangle<Scalar>;
using Bvh = bvh::Bvh<Scalar>;
using BoundingBox = bvh::BoundingBox<Scalar>;
using Sphere = bvh::Sphere<Scalar>;

int main(int argc, char **argv) {
  assert(argc == 2);
  std::string gltf_path(argv[1]);
  std::cerr << "loading " << gltf_path << std::endl;

  tinygltf::Model model;
  tinygltf::TinyGLTF loader;
  std::string err;
  std::string warn;

  bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, gltf_path.c_str());
  if (!warn.empty()) {
    std::cerr << "warning: " << warn;
  }
  if (!err.empty()) {
    std::cerr << "error: " << err;
  }
  if (!ret) {
    std::cerr << "failed to parse glTF" << std::endl;
    return -1;
  }

  assert(model.scenes.size() == 1);
  const auto &scene = model.scenes[0];
  std::map<int, std::set<int>> mesh_nodes;
  std::vector<int> front;
  for (const int nidx : scene.nodes) {
    front.push_back(nidx);
  }
  while (!front.empty()) {
    int nidx = front.back();
    front.pop_back();
    assert(static_cast<size_t>(nidx) < model.nodes.size());
    const auto &node = model.nodes[nidx];
    if (node.mesh != -1) {
      mesh_nodes[node.mesh].insert(nidx);
    }
    for (const int cidx : node.children) {
      front.push_back(cidx);
    }
  }

  std::vector<bvh::Triangle<float>> triangles;
  std::map<int, std::pair<size_t, size_t>> mesh_triangles;
  for (const auto &[midx, nodes] : mesh_nodes) {
    std::ignore = nodes;
    assert(static_cast<size_t>(midx) < model.meshes.size());
    const auto &mesh = model.meshes[midx];
    size_t start = triangles.size();
    for (const auto &prim : mesh.primitives) {
      // TODO: support more than just triangles
      assert(prim.mode == TINYGLTF_MODE_TRIANGLES);

      std::vector<Vec3> vertices;
      {
        const auto pidx = prim.attributes.at("POSITION");
        const tinygltf::Accessor &accessor = model.accessors[pidx];
        const tinygltf::BufferView &bufferView =
            model.bufferViews[accessor.bufferView];
        const tinygltf::Buffer &buffer = model.buffers[bufferView.buffer];
        const int stride = accessor.ByteStride(bufferView);
        assert(accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT &&
               accessor.type == TINYGLTF_TYPE_VEC3);
        const uint8_t *data =
            (const uint8_t *)&buffer
                .data[bufferView.byteOffset + accessor.byteOffset];
        vertices.resize(accessor.count);
        for (size_t i = 0; i < accessor.count; i++) {
          vertices[i] = *reinterpret_cast<const Vec3 *>(&data[i * stride]);
        }
      }

      std::vector<uint32_t> indices;
      {
        const tinygltf::Accessor &accessor = model.accessors[prim.indices];
        const tinygltf::BufferView &bufferView =
            model.bufferViews[accessor.bufferView];
        const tinygltf::Buffer &buffer = model.buffers[bufferView.buffer];
        const int stride = accessor.ByteStride(bufferView);
        assert(
            (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE ||
             accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT ||
             accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) &&
            accessor.type == TINYGLTF_TYPE_SCALAR);

        indices.resize(accessor.count);
        const uint8_t *data =
            (const uint8_t *)&buffer
                .data[bufferView.byteOffset + accessor.byteOffset];
        for (size_t i = 0; i < accessor.count; i++) {
          switch (accessor.componentType) {
          case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE: {
            const uint8_t *idx = (uint8_t *)&data[i * stride];
            indices[i] = static_cast<uint32_t>(idx[0]);
            break;
          }
          case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT: {
            const uint16_t *idx = (uint16_t *)&data[i * stride];
            indices[i] = static_cast<uint32_t>(idx[0]);
            break;
          }
          case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT: {
            const uint32_t *idx = (uint32_t *)&data[i * stride];
            indices[i] = static_cast<uint32_t>(idx[0]);
            break;
          }
          default:
            assert(false);
          }
        }
      }
      assert(indices.size() % 3 == 0);
      for (size_t i = 0; i < indices.size(); i += 3) {
        triangles.push_back(Triangle(vertices[indices[i]],
                                     vertices[indices[i + 1]],
                                     vertices[indices[i + 2]]));
      }
    }
    size_t end = triangles.size() - 1;
    mesh_triangles[midx] = {start, end};
  }

  std::cerr << "found " << triangles.size() << " triangles. constructing BVH"
            << std::endl;
  Bvh bvh;
  auto [bboxes, centers] = bvh::compute_bounding_boxes_and_centers(
      triangles.data(), triangles.size());
  auto global_bbox =
      bvh::compute_bounding_boxes_union(bboxes.get(), triangles.size());
  bvh::SweepSahBuilder<Bvh> builder(bvh);
  builder.build(global_bbox, bboxes.get(), centers.get(), triangles.size());

  Scalar radius = bvh::length(global_bbox.diagonal());
  Sphere bsphere(global_bbox.center(), radius);
  std::cerr << "bounding sphere at (" << bsphere.origin[0] << ", "
            << bsphere.origin[1] << ", " << bsphere.origin[2] << ") w/ radius "
            << bsphere.radius << std::endl;

  // TODO: Sun is fixed at high noon right now. Make this configurable.
  Vec3 sun_center = bsphere.origin + Vec3(0, 0, bsphere.radius);
  std::cerr << "sun disk at (" << sun_center[0] << ", " << sun_center[1] << ", "
            << sun_center[2] << ")" << std::endl;

  // First pass to accumulate light on each mesh
  for (const auto &[midx, interval] : mesh_triangles) {
    const auto &mesh = model.meshes[midx];
    const auto [start, end] = interval;
    std::cerr << "launching " << end - start << " rays from " << mesh.name
              << std::endl;

    size_t hits = 0;
    size_t self_hits = 0;
    size_t misses = 0;
    for (size_t i = start; i <= end; i++) {
      const auto &t = triangles[i];
      const Vec3 center = t.center();

      // TODO: cast to closest point on disk instead of center
      bvh::Ray<Scalar> ray(center, center - sun_center, 0.0, bsphere.radius);
      bvh::ClosestPrimitiveIntersector<Bvh, Triangle> primitive_intersector(
          bvh, triangles.data());
      bvh::SingleRayTraverser<Bvh> traverser(bvh);

      if (auto hit = traverser.traverse(ray, primitive_intersector)) {
        auto tidx = hit->primitive_index;
        if (tidx >= start && tidx <= end) {
          // Hit mesh we originated from
          self_hits++;
        } else {
          // Hit another mesh
          hits++;
        }
      } else {
        // Miss! Saw the sun
        misses++;
      }
    }
    std::cerr << "\thits: " << hits << std::endl
              << "\tself_hits: " << self_hits << std::endl
              << "\tmisses: " << misses << std::endl;
  }

  // TODO: Second pass to scatter some portion of light form each surface

  // TODO: output copy of gltf from user w/ radiance annotations

  return 0;
}
