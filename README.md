# Tracer

`tracer` identifies objects in a scene and does path tracing to measure how
exposed they are to sunlight.

This is more-or-less an implementation of ["A reverse ray-tracing method for modelling the net radiative flux in leaf-resolving plant canopy simulations"](https://doi.org/10.1016/j.ecolmodel.2017.11.022)
by Brian N. Bailey.

## Dependencies

You'll need `meson` and `ninja` to build. You can get them from `pip` or your
favorite packager manager.
```
pip install meson ninja
```

- (Required) [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
- (Required) [GLM](https://github.com/g-truc/glm) 
- (Optional) [OpenMP](https://www.openmp.org/)

On Debian you can install these dependencies with:
```
apt-get install libomp-dev libglm-dev pkg-config
```

## Building
To build, run:
```
git clone <this repo>
cd tracer
git submodule update --init
meson setup build
ninja -C build
```

## Supported Input Formats
- (DONE) GLTF
- (TODO) OBJ
- (TODO) PLY

## Usage
`tracer` takes the path to the scene it should measure as its sole command line
argument.
