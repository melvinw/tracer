# Tracer

`tracer` identifies objects in a scene and does path tracing to measure how
exposed they are to sunlight.

## Dependencies

You'll need `meson` and `ninja` to build. You can get them from `pip` or your
favorite packager manager.
```
pip install meson ninja
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
