#ifndef GLTF_LOADER_H_
#define GLTF_LOADER_H_

#include <utility>
#include <vector>

#include "types.h"

using MeshList = std::vector<Mesh>;
using MeshInstances = std::vector<MeshInstance>;
std::pair<MeshList, MeshInstances> LoadMeshesFromGLTF(const std::string &path);

#endif // GLTF_LOADER_H_
