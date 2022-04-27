#include "types.h"

using MeshList = std::vector<Mesh>;
using MeshInstances = std::vector<MeshInstance>;
std::pair<MeshList, MeshInstances> LoadMeshesFromOBJ(const std::string &path);