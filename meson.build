project('tracer', 'cpp')

inc_dirs = [
	include_directories(
    'thirdparty/CLI11/include/',
		'thirdparty/bvh/include',
		'thirdparty/thread-pool',
		'thirdparty/tinygltf',
    'thirdparty/tinyobjloader'
	),
	include_directories(
		'include',
		is_system: false
	),
]

cxx_flags = [
  '-std=c++17',
  '-Wall',
  '-Wextra',
  '-Werror',
  '-Wno-error=unknown-pragmas'
]

omp_dep = dependency('openmp', required: false)
glm_dep = dependency('glm', required: true)
threads_dep = dependency('threads', required: true)

tinygltf = build_target(
  'tinygltf',
  'thirdparty/tinygltf/tiny_gltf.cc',
  target_type : 'library',
  cpp_args: cxx_flags
)

tinyobjloader = build_target(
  'tinyobjloader',
  'thirdparty/tinyobjloader/tiny_obj_loader.cc',
  target_type : 'library',
  cpp_args: cxx_flags
)


executable(
  'tracer',
  ['src/main.cc', 'src/gltf_loader.cc', 'src/gpl/solar_position.cc', 'src/obj_loader.cc'],
  include_directories: inc_dirs,
  link_with: [tinygltf,tinyobjloader],
  dependencies: [glm_dep, omp_dep, threads_dep],
  cpp_args: cxx_flags
)
