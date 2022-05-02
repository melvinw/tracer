#ifndef TYPES_H_
#define TYPES_H_

#include <cinttypes>
#include <optional>
#include <string>
#include <vector>

#include <glm/glm.hpp>

struct Vertex {
  glm::vec4 position;
  glm::vec4 normal;
};

struct Mesh {
  std::vector<Vertex> vertices;
  std::vector<size_t> indices;
  glm::vec3 min;
  glm::vec3 max;
};

struct MeshInstance {
  std::string name;
  size_t mesh_index;
  std::optional<glm::mat4> model_matrix;
};

using MeshList = std::vector<Mesh>;
using MeshInstances = std::vector<MeshInstance>;

#endif  // TYPES_H_
