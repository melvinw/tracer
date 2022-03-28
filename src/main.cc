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

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/transform.hpp>
#include <tiny_gltf.h>

using Scalar = double;
using Vec3 = bvh::Vector3<Scalar>;
using Triangle = bvh::Triangle<Scalar>;
using Bvh = bvh::Bvh<Scalar>;
using BoundingBox = bvh::BoundingBox<Scalar>;
using Sphere = bvh::Sphere<Scalar>;

struct Vertex {
  glm::vec4 position;
  glm::vec4 normal;
};

struct Mesh {
  std::vector<Vertex> vertices;
  std::vector<uint32_t> indices;
};

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

  std::vector<Mesh> meshes(model.meshes.size());
  for (const auto &[midx, nodes] : mesh_nodes) {
    std::ignore = nodes;
    assert(static_cast<size_t>(midx) < model.meshes.size());
    const auto &mesh = model.meshes[midx];
    for (const auto &prim : mesh.primitives) {
      // TODO: support more than just triangles
      assert(prim.mode == TINYGLTF_MODE_TRIANGLES);

      auto &vertices = meshes[midx].vertices;
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
          vertices[i].position = glm::vec4(
              reinterpret_cast<const glm::vec3 *>(&data[i * stride])[0], 1);
        }
      }

      if (prim.attributes.count("NORMAL")) {
        const auto nidx = prim.attributes.at("NORMAL");
        const tinygltf::Accessor &accessor = model.accessors[nidx];
        const tinygltf::BufferView &bufferView =
            model.bufferViews[accessor.bufferView];
        const tinygltf::Buffer &buffer = model.buffers[bufferView.buffer];
        const int stride = accessor.ByteStride(bufferView);
        assert(accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT &&
               accessor.type == TINYGLTF_TYPE_VEC3);
        const uint8_t *data =
            (const uint8_t *)&buffer
                .data[bufferView.byteOffset + accessor.byteOffset];
        assert(accessor.count == vertices.size());
        for (size_t i = 0; i < accessor.count; i++) {
          vertices[i].normal = glm::vec4(
              reinterpret_cast<const glm::vec3 *>(&data[i * stride])[0], 0);
        }
      }

      auto &indices = meshes[midx].indices;
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
    }
  }

  std::vector<bvh::Triangle<Scalar>> triangles;
  std::map<int, std::pair<size_t, size_t>> node_triangles;
  for (const auto &[midx, nodes] : mesh_nodes) {
    const auto &mesh = meshes[midx];
    for (const auto nidx : nodes) {
      const auto &node = model.nodes[nidx];
      size_t start = triangles.size();

      glm::vec3 translation =
          node.translation.size() == 3
              ? glm::vec3(node.translation[0], node.translation[1],
                          node.translation[2])
              : glm::vec3(0.0f);
      glm::vec4 rotation = node.rotation.size() == 4
                               ? glm::vec4(node.rotation[0], node.rotation[1],
                                           node.rotation[2], node.rotation[3])
                               : glm::vec4(0.0f);
      glm::vec3 scale =
          node.scale.size() == 3
              ? glm::vec3(node.scale[0], node.scale[1], node.scale[2])
              : glm::vec3(1.0f);
      auto trs = glm::mat4(1.0f);
      if (glm::any(glm::notEqual(translation, glm::vec3(0.0f)))) {
        trs = glm::translate(trs, translation);
      }
      if (glm::any(glm::notEqual(rotation, glm::vec4(0.0f)))) {
        trs *= glm::toMat4(
            glm::quat(rotation.w, rotation.x, rotation.y, rotation.z));
      }
      if (glm::any(glm::notEqual(scale, glm::vec3(1.0f)))) {
        trs = glm::scale(trs, scale);
      }

      for (size_t i = 0; i < mesh.indices.size(); i += 3) {
        glm::vec4 v0 = trs * mesh.vertices[mesh.indices[i]].position;
        glm::vec4 v1 = trs * mesh.vertices[mesh.indices[i + 1]].position;
        glm::vec4 v2 = trs * mesh.vertices[mesh.indices[i + 2]].position;
        triangles.push_back(Triangle(Vec3(v0.x, v0.y, v0.z),
                                     Vec3(v1.x, v1.y, v1.z),
                                     Vec3(v2.x, v2.y, v2.z)));
      }
      size_t end = triangles.size() - 1;
      node_triangles[nidx] = {start, end};
    }
  }

  std::cerr << "found " << model.nodes.size() << " models consisting of "
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
  for (const auto &[nidx, interval] : node_triangles) {
    const auto &node = model.nodes[nidx];
    const auto [start, end] = interval;
    std::cerr << "launching " << end - start + 1 << " rays from " << node.name
              << "(" << start << ", " << end << ")" << std::endl;

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
