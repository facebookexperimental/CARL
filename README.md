# CARL

CARL -- an initialism of Classical Action Recognition Library -- is a 
C++ library for generalized action recognition, primarily intended to 
enable power "gesture recognition" in XR applications.

## Getting Started

- [Standard Repo Set Up](#standard-repo-set-up)
- [Building for WebAssembly with Emscripten](#building-for-webassembly-with-emscripten)
- [Using CARL](#using-carl)

### Standard Repo Set Up

CARL uses CMake for build system generation and Git submodules to
provide dependencies. Substituting other build systems and dependency
approaches, such as Buck, is fairly straightforward but is not 
discussed in this document.

1. Clone the repo:
   ```
   git clone https://github.com/facebookexperimental/CARL.git
   ```
1. Navigate into the repo:
   ```
   cd CARL
   ```
1. Install dependencies: 
   ```
   git submodule update --init --recursive
   ```
1. Create and navigate to the build directory:
   ```
   cd ..
   mkdir CARL_build
   cd CARL_build
   ```
1. Invoke CMake to generate your preferred build system. For example, on
   Windows with VS2019 installed, default arguments will generate a Visual
   Studio Solution called `carl.sln`.
   ```
   cmake ../CARL
   ```

### Building for Android on with the NDK

These instructions are known to work on Linux. Note that this only works on 
a full real Linux install; the NDK CMake integration does not currently seem
to work well on Windows Subsystem for Linux, unfortunately.

1. Download and unzip the [NDK](https://developer.android.com/ndk/downloads).
1. Set up a CARL repository as described [above](@standard-repo-set-up),
   omitting the build system creation.
1. Invoke CMake using 
   [CMake's built-in Android integration](https://cmake.org/cmake/help/latest/manual/cmake-toolchains.7.html#cross-compiling-for-android).
   If you need to choose Android platform settings other than the 
   defaults, they should likely be specified in this step.
   ```
   cmake ../CARL -D CMAKE_SYSTEM_NAME=Android -D CMAKE_BUILD_TYPE=Release -D CMAKE_ANDROID_ARCH_ABI="arm64-v8a" -D CMAKE_ANDROID_NDK=[NDK root directory]
   ```
   **Note:** use of the the `CMAKE_BUILD_TYPE` and `CMAKE_ANDROID_ARCH_ABI` 
   variables as shown above is recommended when building CARL for use as a 
   Unity plugin on the Meta Quest 3; in particular, neglecting to specify 
   an ABI matching that of your app can cause Unity to fail to load the 
   library at runtime, which can be tricky to diagnose. Strictly speaking,
   however, specifying these variables is not required and they can be 
   omitted in situations where their default values are appropriate.
1. Build the `carl_shared_library` target to produce a `.so` which can be 
   used in Android applications.
   ```
   make carl_shared_library
   ```

### Building for WebAssembly with Emscripten

1. Set up Emscripten using the 
   [recommended emsdk instructions](https://emscripten.org/docs/getting_started/downloads.html#installation-instructions-using-the-emsdk-recommended).
   On Windows, it is highly recommended to do this using the
   [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) 
   to set up the build in a Unix-like environment.
1. Set up a CARL repository as described [above](#standard-repo-set-up),
   omitting the build system creation.
1. Invoke CMake using the Emscripten toolchain to generate a WebAssembly
   build system.
   ```
   cmake ../CARL -D CMAKE_TOOLCHAIN_FILE=[TODO: Insert directory path here]/Emscripten.cmake
   ```
1. Build the `carl_static_library` target to produce the WebAssembly 
   artifacts. Depending on your Shell's permissions, you may need to 
   add `sudo` as follows:
   ```
   emmake sudo make carl_static_library
   ```

### Using CARL

TODO: Information coming soon.

## Dependencies and Licensing

CARL itself is licensed under the MIT license. CARL has three source 
dependencies:

- [arcana.cpp](https://github.com/microsoft/arcana.cpp): Arcana is a 
  lightweight C++ utility library published and maintained by Microsoft.
  It is made available under the MIT License.
- [GSL](https://github.com/microsoft/gsl): The Guidelines Support Library
  is another lightweight C++ utility library published and maintained by 
  Microsoft. It is also available under the MIT License.
- [Eigen](https://gitlab.com/libeigen/eigen): Eigen is an advanced 
  mathematics library providing various linear algebra and 3D arithmetic
  features. It is available under several licenses, but CARL specifically
  depends on the MPL2 subset and so uses Eigen under the Mozilla Public 
  License 2.0.

## Contributing

Development of CARL happens in the open on GitHub, and we are grateful to 
the community for contributing bugfixes and improvements. Read below to 
learn how you can take part in improving CARL.

- [Code of Conduct](CODE_OF_CONDUCT.md)
- [Contributing Guide](CONTRIBUTING.md)
