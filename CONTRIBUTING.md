# Contributing

## Code Style

This code in this repo mostly conforms to the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html). Try to adhere to it when submitting changes.

All changes should be run through `clang-format` version 11 before review, to minimize style-related quibbling.

For example, you can format the entire project by running the following.
```
clang-format-11 -i include/**/*.h src/**/*.cc
```
You can find some documentation about `clang-format-11` [here](https://releases.llvm.org/11.0.0/tools/clang/docs/ClangFormat.html).

If your package manager of choice doesn't have `clang-format-11` you can find prebuilt binaries for the following distributions at the links below:
- [macos](https://github.com/llvm/llvm-project/releases/download/llvmorg-11.0.0/clang+llvm-11.0.0-x86_64-apple-darwin.tar.xz)
- [Ubuntu 16.04](https://github.com/llvm/llvm-project/releases/download/llvmorg-11.0.0/clang+llvm-11.0.0-x86_64-linux-gnu-ubuntu-16.04.tar.xz)
- [Ubuntu 20.04](https://github.com/llvm/llvm-project/releases/download/llvmorg-11.0.0/clang+llvm-11.0.0-x86_64-linux-gnu-ubuntu-20.04.tar.xz)
- [FreeBSD 11](https://github.com/llvm/llvm-project/releases/download/llvmorg-11.0.0/clang+llvm-11.0.0-amd64-unknown-freebsd11.tar.xz)

For other platforms/OSes, see https://releases.llvm.org/download.html.

## Dependencies

Ideally the project will depend on as few things as possible (especially at runtime), but sometimes it makes sense to fold in a library to do the heavy lifting for a feature. These dependencies are kept in the `thirdparty` folder. How they are incorporated into the build system will vary on a case-by-case basis. See `meson.build` for hints about how other dependencies have been integrated.

Make sure any dependencies you add have licenses that are compatible with the license in `LICENSE`.

## Dependencies Distributed with Git

If the dependency you intend to add is tracked in a git repo, add it as a submodule using the following command.
```
git submodule add <repo_url> thirdparty/<dependency_name>
```

## Dependencies Distributed as Archives

If the dependency you intend to add is distributed as an archive download and unpack it directly into `thirdpary/<dependency_name>` and track it with `git add thirdparty/<dependency_name>`.
