# CARL — Classical Action Recognition Library

CARL is a C++ library for real-time action recognition in XR applications. Instead of relying on machine learning models or large datasets, CARL uses classical techniques to recognize actions defined purely by example — just perform the gesture a few times and CARL learns to recognize it.

## Key Features

- **Define by example** — no training data, no ML expertise; perform an action a few times to create a definition
- **Real-time recognition** — continuous scoring with sub-frame latency on mobile XR hardware
- **12 built-in action types** — hand poses, hand gestures, controller gestures, wrist trajectories, and more
- **Custom action types** — extend CARL with your own descriptor logic
- **Counterexamples** — add negative examples to disambiguate similar actions
- **Sensitivity tuning** — per-definition and per-recognizer sensitivity controls
- **Auto-trimming** — automatically determine optimal example boundaries
- **Cross-platform** — Windows, Android (Quest), WebAssembly; Unity integration included
- **No external runtime dependencies** — pure C++ with a C API for FFI

## Quick Start (C++)

```cpp
#include <carl/Carl.h>
#include <carl/utilities/FileSerialization.h>
#include <iostream>

using namespace carl;
using namespace carl::action;

int main() {
    // 1. Create a session to receive live input
    Session session{/* singleThreaded */ true};

    // 2. Record some input samples (in practice, these come from XR hardware)
    InProgressRecording inProgress{/* maxSeconds */ 10};
    for (int i = 0; i < 100; ++i) {
        InputSample sample{};
        sample.Timestamp = i * 0.016; // ~60 fps
        // Populate sample.LeftHandJointPoses, sample.HmdPose, etc.
        inProgress.addSample(std::move(sample));
    }
    Recording recording{std::move(inProgress)};

    // 3. Create an example marking when the action occurred
    double startTime = 0.5;
    double endTime = 1.2;
    Example example{recording, startTime, endTime};

    // 4. Build a definition with at least one example
    Definition definition{ActionType::RightHandGesture};
    definition.addExample(std::move(example));

    // 5. Create a recognizer and listen for recognition events
    Recognizer recognizer{session, definition};

    auto ticket = recognizer.whenRecognitionChangedSignal.addHandler(
        [](bool recognized) {
            std::cout << (recognized ? "Action recognized!" : "Action ended.") << std::endl;
        });

    // 6. Feed live input and tick callbacks
    arcana::cancellation token{};
    for (/* each frame */) {
        InputSample liveSample{};
        liveSample.Timestamp = /* current time in seconds */;
        // Populate liveSample fields...
        session.addInput(liveSample);
        session.tickCallbacks(token);
    }

    return 0;
}
```

## Documentation

| Document | Description |
|----------|-------------|
| [Concepts](Documentation/Concepts.md) | Core concepts: InputSample, Recording, Example, Definition, Session, Recognizer |
| [Action Types](Documentation/ActionTypes.md) | All 12 action types explained with selection guidance |
| [Defining Actions](Documentation/Tutorials/DefiningActions.md) | Tutorial: recording input and building definitions |
| [Recognizing Actions](Documentation/Tutorials/RecognizingActions.md) | Tutorial: creating sessions, recognizers, and handling events |
| [Tuning and Debugging](Documentation/Tutorials/TuningAndDebugging.md) | Sensitivity, counterexamples, auto-trimming, and troubleshooting |
| [C++ API Reference](Documentation/CppApiReference.md) | Complete C++ API reference |
| [C API Reference](Documentation/CApiReference.md) | Complete C API reference (41 functions) |
| [Serialization](Documentation/Serialization.md) | Saving and loading definitions and recordings |
| [FAQ](Documentation/FAQ.md) | Common questions and migration notes |
| [Unity Integration](Targets/unity/README.md) | Unity package setup and usage |

## Building CARL

CARL uses CMake for build system generation and Git submodules for dependencies.

### Standard Setup

1. Clone the repo:
   ```
   git clone https://github.com/facebookexperimental/CARL.git
   ```
2. Navigate into the repo:
   ```
   cd CARL
   ```
3. Install dependencies:
   ```
   git submodule update --init --recursive
   ```
4. Create and navigate to the build directory:
   ```
   cd ..
   mkdir CARL_build
   cd CARL_build
   ```
5. Invoke CMake to generate your preferred build system. For example, on
   Windows with VS2019 installed, default arguments will generate a Visual
   Studio Solution called `carl.sln`.
   ```
   cmake ../CARL
   ```

### Building for Android with the NDK

#### Building With Ninja on Windows

By far the easiest way to build for Android on Windows is using [Ninja](https://ninja-build.org/);
simply have a Ninja binary in a directory accessible from your PATH, and CMake will be able to find
it and configure a build from a command like the following:

```
cmake -G Ninja ../CARL -D CMAKE_SYSTEM_NAME=Android -D CMAKE_ANDROID_ARCH_ABI="arm64-v8a" -D CMAKE_ANDROID_NDK=C:/Microsoft/AndroidNDK/android-ndk-r23c -D CMAKE_ANDROID_STL_TYPE=c++_static
```

This command uses an NDK directory of `C:/Microsoft/AndroidNDK/android-ndk-r23c`, which is where
Visual Studio may put an NDK [if you ask it to install one for you](https://learn.microsoft.com/en-us/windows/android/native-android#use-c-or-c-for-android-game-development).
If you have the NDK installed elsewhere, simply use that directory path in the command above instead.

Once the build is configured, compiling is as simple as invoking

```
ninja carl_shared_library
```

#### Building With Make on Linux

Note that this only works on a full real Linux install; the NDK CMake integration does not
currently seem to work well on Windows Subsystem for Linux, unfortunately.

1. Download and unzip the [NDK](https://developer.android.com/ndk/downloads).
2. Set up a CARL repository as described [above](#standard-setup),
   omitting the build system creation.
3. Invoke CMake using
   [CMake's built-in Android integration](https://cmake.org/cmake/help/latest/manual/cmake-toolchains.7.html#cross-compiling-for-android).
   If you need to choose Android platform settings other than the
   defaults, they should likely be specified in this step.
   ```
   cmake ../CARL -D CMAKE_SYSTEM_NAME=Android -D CMAKE_BUILD_TYPE=Release -D CMAKE_ANDROID_ARCH_ABI="arm64-v8a" -D CMAKE_ANDROID_NDK=[NDK root directory]
   ```
   **Note:** use of the `CMAKE_BUILD_TYPE` and `CMAKE_ANDROID_ARCH_ABI`
   variables as shown above is recommended when building CARL for use as a
   Unity plugin on the Meta Quest 3; in particular, neglecting to specify
   an ABI matching that of your app can cause Unity to fail to load the
   library at runtime, which can be tricky to diagnose. Strictly speaking,
   however, specifying these variables is not required and they can be
   omitted in situations where their default values are appropriate.
4. Build the `carl_shared_library` target to produce a `.so` which can be
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
2. Set up a CARL repository as described [above](#standard-setup),
   omitting the build system creation.
3. Invoke CMake using the Emscripten toolchain to generate a WebAssembly
   build system.
   ```
   cmake ../CARL -D CMAKE_TOOLCHAIN_FILE=[emsdk directory]/Emscripten.cmake
   ```
4. Build the `carl_static_library` target to produce the WebAssembly
   artifacts. Depending on your shell's permissions, you may need to
   add `sudo` as follows:
   ```
   emmake make carl_static_library
   ```

## Coordinate System

CARL uses a **right-handed Y-up** coordinate system, consistent with the [OpenXR](https://registry.khronos.org/OpenXR/specs/1.0/html/xrspec.html) convention. All 3D poses in `InputSample` must be expressed in this coordinate system.

If you are using CARL in Unity (which uses left-handed Y-up), the `CarlXRInputProvider` component handles the conversion automatically. For manual conversion, negate the Z component of positions and the Z and W components of quaternions. See the [Unity Integration](Targets/unity/README.md) docs for details.

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
