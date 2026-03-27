# CARL Unity Integration

Unity package for the **Classical Action Recognition Library (CARL)** — a gesture recognition system for XR applications that uses Dynamic Time Warping (DTW) to match gestures defined by example.

## Requirements

- Unity 2021.3+
- Native CARL plugin binary (`carl.dll` for Windows, `libcarl.so` for Android/Quest)
- Optional: [Unity XR Hands](https://docs.unity3d.com/Packages/com.unity.xr.hands@1.3) package (≥1.3.0) for hand tracking input

## Installation

1. Build the native CARL library for your target platform (see [Building Native Plugins](#building-native-plugins)).
2. Add this package to your Unity project via the Package Manager:
   - **Local path**: Point to `Targets/unity/` in the CARL repository.
   - **Git URL**: Use the repository URL with `?path=Targets/unity`.
3. Place the native plugin binaries in the appropriate `Plugins/` directory.

## Quick Start

1. Add a **CarlSessionManager** component to a GameObject in your scene.
2. Add a **CarlXRInputProvider** component (same or different GameObject).
3. Create a `CarlDefinitionAsset` (right-click in Project → Create → CARL → Definition Asset) and import a `.carl` definition file.
4. Add a **CarlGestureRecognizer** component and assign your definition asset.
5. Wire up the `OnGestureDetected` / `OnGestureLost` Unity Events.

## Architecture

```
┌─────────────────────────────────────┐
│         Unity Components            │
│  CarlSessionManager                 │
│  CarlXRInputProvider                │
│  CarlGestureRecognizer              │
│  CarlGestureRecorder                │
├─────────────────────────────────────┤
│       Managed Wrapper Layer         │
│  CarlSession    CarlDefinition      │
│  CarlRecognizer CarlRecording       │
│  CarlExample    CarlRecordingInsp.  │
├─────────────────────────────────────┤
│        Native P/Invoke Layer        │
│  CarlNative  CarlInteropTypes      │
│  CarlBytes                          │
├─────────────────────────────────────┤
│     Native CARL Library (C API)     │
│  carl.dll / libcarl.so              │
└─────────────────────────────────────┘
```

## Coordinate System

CARL uses the **OpenXR right-handed Y-up** coordinate system. Unity uses **left-handed Y-up**. The `CarlXRInputProvider` handles conversion automatically. If you feed input samples manually, use `CarlCoordinateConversion` to convert poses.

## Threading Model

- `CarlSession` runs gesture processing on a background thread by default.
- `CarlSessionManager.Update()` calls `TickCallbacks()` each frame to dispatch results to the main thread.
- `carl_getCurrentScore()` is thread-safe and can be polled at any time.

## Building Native Plugins

### Windows (x64)
```bash
mkdir build && cd build
cmake ../CARL
cmake --build . --config Release
# Copy carl.dll to Targets/unity/Runtime/Plugins/x86_64/
```

### Android (arm64-v8a, for Quest)
```bash
mkdir build-android && cd build-android
cmake -G Ninja ../CARL \
  -DCMAKE_SYSTEM_NAME=Android \
  -DCMAKE_ANDROID_ARCH_ABI="arm64-v8a" \
  -DCMAKE_ANDROID_NDK=/path/to/ndk \
  -DCMAKE_ANDROID_STL_TYPE=c++_static \
  -DCMAKE_BUILD_TYPE=Release
ninja carl
# Copy libcarl.so to Targets/unity/Runtime/Plugins/Android/arm64-v8a/
```

## License

MIT — see [LICENSE.md](LICENSE.md).
