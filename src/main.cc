#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include <bvh/bvh.hpp>
#include <bvh/sweep_sah_builder.hpp>
#include <bvh/triangle.hpp>
#include <bvh/vector.hpp>
#include <tiny_gltf.h>

using Vec3 = bvh::Vector3<float>;
using Triangle = bvh::Triangle<float>;
using Bvh = bvh::Bvh<float>;

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
  for (const int nidx : scene.nodes) {
    assert(static_cast<size_t>(nidx) < model.nodes.size());
    const auto &node = model.nodes[nidx];
    if (node.mesh != -1) {
      mesh_nodes[node.mesh].insert(nidx);
    }
    // TODO: allow more than one level of nesting
    for (const int cidx : node.children) {
      assert(static_cast<size_t>(cidx) < model.nodes.size());
      const auto &child = model.nodes[cidx];
      if (child.mesh != -1) {
        mesh_nodes[child.mesh].insert(cidx);
      }
    }
  }

  std::vector<bvh::Triangle<float>> triangles;
  for (const auto &[midx, nodes] : mesh_nodes) {
    std::ignore = nodes;
    assert(static_cast<size_t>(midx) < model.meshes.size());
    const auto &mesh = model.meshes[midx];
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

  return 0;
}
