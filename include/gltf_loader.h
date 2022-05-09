#ifndef GLTF_LOADER_H_
#define GLTF_LOADER_H_

#include <utility>
#include <vector>

#include "types.h"

std::pair<MeshList, MeshInstances> LoadMeshesFromGLTF(const std::string &path);

#endif  // GLTF_LOADER_H_
