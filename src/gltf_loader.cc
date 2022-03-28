#include "gltf_loader.h"

#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/transform.hpp>

#include <tiny_gltf.h>

std::pair<MeshList, MeshInstances> LoadMeshesFromGLTF(const std::string &path) {
  tinygltf::Model model;
  tinygltf::TinyGLTF loader;
  std::string err;
  std::string warn;

  bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, path.c_str());
  if (!warn.empty()) {
    std::cerr << "warning: " << warn;
  }
  if (!err.empty()) {
    std::cerr << "error: " << err;
  }
  if (!ret) {
    std::cerr << "failed to parse glTF" << std::endl;
    assert(false);
  }

  assert(model.scenes.size() == 1);
  const auto &scene = model.scenes[0];
  std::map<size_t, std::set<size_t>> mesh_nodes;
  std::vector<size_t> front;
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
    assert(midx < model.meshes.size());
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

  std::vector<MeshInstance> instances;
  for (const auto &[midx, nodes] : mesh_nodes) {
    // TODO: accumulate transforms through ancestor graph
    for (const auto nidx : nodes) {
      const auto &node = model.nodes[nidx];
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
      instances.push_back(MeshInstance{
        name : node.name,
        mesh_index : midx,
        model_matrix : trs
      });
    }
  }

  return std::make_pair(meshes, instances);
}
