project('tracer', 'cpp')

inc_dirs = include_directories(
  'thirdparty/bvh/include',
  'thirdparty/tinygltf'
)
cxx_flags = [
  '-std=c++17',
  '-Wall',
  '-Wextra',
  '-Werror',
  '-Wno-error=unknown-pragmas'
]

tinygltf = build_target(
  'tinygltf',
  'thirdparty/tinygltf/tiny_gltf.cc',
  target_type : 'library',
  cpp_args: cxx_flags
)

executable(
  'tracer',
  'src/main.cc',
  include_directories: inc_dirs,
  link_with: tinygltf,
  cpp_args: cxx_flags
)
