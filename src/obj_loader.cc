#include "obj_loader.h"
#include "tiny_obj_loader.h"

#include <iostream>
#include <set>
std::pair<MeshList, MeshInstances> LoadMeshesFromOBJ(const std::string &path) {
  std::cout << "Test" << std::endl;
  bool ret = false;
  tinyobj::attrib_t attribute;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;
  std::string warn;
  std::string err;
  ret = tinyobj::LoadObj(&attribute, &shapes, &materials, &warn, &err,
                         path.c_str(), NULL, true,
                         false);  // TODO: check what last parameter means
  if (!ret) {
    std::cerr << warn << " " << err << std::endl;
  }

  MeshList meshes(shapes.size());
  MeshInstances meshInstances(shapes.size());
  for (size_t shape_index = 0; shape_index < shapes.size(); shape_index++) {
    auto shape = shapes[shape_index];  // Unnecessary copying?
    std::set<std::pair<int, int>> unique_indices;
    std::vector<size_t> indices;
    for (auto indice : shape.mesh.indices) {
      unique_indices.insert(
          std::make_pair(indice.vertex_index, indice.normal_index));
      indices.push_back(indice.vertex_index);
    }
    std::cout << indices.size() << std::endl;
    std::vector<Vertex> vertices;
    for (auto indice : unique_indices) {
      Vertex v;
      v.position = glm::vec4(attribute.vertices[indice.first * 3],
                             attribute.vertices[indice.first * 3 + 1],
                             attribute.vertices[indice.first * 3 + 2],
                             1);  // TODO: last one taken as 1
      v.normal = glm::vec4(attribute.normals[indice.second * 3],
                           attribute.normals[indice.second * 3 + 1],
                           attribute.normals[indice.second * 3 + 2],
                           1);  // TODO: handle the case when there are no
                                // normals and we have to generate them
      vertices.push_back(v);
    }

    meshes[shape_index].vertices = vertices;
    meshes[shape_index].indices = indices;
    meshInstances[shape_index].name = shape.name;
    meshInstances[shape_index].mesh_index = shape_index;
    meshInstances[shape_index].model_matrix = glm::mat4(1.0f);
  }
  return make_pair(meshes, meshInstances);
}