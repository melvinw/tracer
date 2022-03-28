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
};

struct MeshInstance {
  std::string name;
  size_t mesh_index;
  std::optional<glm::mat4> model_matrix;
};

#endif // TYPES_H_
