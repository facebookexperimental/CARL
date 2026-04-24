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

CARL uses the **OpenXR right-handed Y-up** coordinate system. Unity uses **left-handed Y-up**. The `CarlXRInputProvider` handles conversion automatically. If you feed input samples manually, use `CarlCoordinateConversion` to convert poses:

```csharp
// Unity → CARL
CarlOptionalTransform carlPose = CarlCoordinateConversion.ToCarlTransform(unityPose);

// CARL → Unity
Pose unityPose = CarlCoordinateConversion.ToUnityPose(carlTransform);
```

Conversion negates the Z component of positions and the Z and W components of quaternions.

## Threading Model

- `CarlSession` runs gesture processing on a background thread by default.
- `CarlSessionManager.Update()` calls `TickCallbacks()` each frame to dispatch results to the main thread.
- `carl_getCurrentScore()` is thread-safe and can be polled at any time.

## Components

### CarlSessionManager

Singleton MonoBehaviour that owns the `CarlSession`. Add one to your scene.

| Inspector Field | Default | Description |
|----------------|---------|-------------|
| `enableLogging` | `true` | Routes CARL diagnostics to `Debug.Log` |
| `singleThreaded` | `false` | Debug mode: all processing inline on main thread |

Access the session via `CarlSessionManager.Instance.Session`.

### CarlXRInputProvider

Reads XR hardware each frame and feeds InputSamples to the session.

| Inspector Field | Default | Description |
|----------------|---------|-------------|
| `captureHands` | `true` | Capture hand joint data from XRHandSubsystem |
| `captureControllers` | `true` | Capture controller buttons, triggers, thumbstick |
| `captureHmd` | `true` | Capture HMD head pose |
| `palmWristSwapWorkaround` | `true` | Workaround for Unity XRHandSubsystem Palm/Wrist swap bug |

### CarlGestureRecognizer

High-level gesture detection with threshold hysteresis.

| Inspector Field | Default | Description |
|----------------|---------|-------------|
| `definitionAsset` | — | `CarlDefinitionAsset` ScriptableObject |
| `definitionFilePath` | — | Fallback: load definition from file path |
| `activationThreshold` | `0.8` | Score must exceed this to fire `OnGestureDetected` |
| `deactivationThreshold` | `0.6` | Score must fall below this to fire `OnGestureLost` |
| `sensitivity` | `-1` | Override sensitivity (-1 = use definition default) |

**Events:**

| Event | Description |
|-------|-------------|
| `OnScoreChanged(float)` | Fires every frame with the current score |
| `OnGestureDetected` | Fires when score crosses above `activationThreshold` |
| `OnGestureLost` | Fires when score drops below `deactivationThreshold` |

**Properties:**

| Property | Type | Description |
|----------|------|-------------|
| `IsRecognized` | `bool` | Whether the gesture is currently detected |
| `LastScore` | `float` | Most recent score value |

The hysteresis (separate activation/deactivation thresholds) prevents flickering when the score hovers near the threshold.

### CarlGestureRecorder

Records gestures at runtime to create new training examples.

| Event | Description |
|-------|-------------|
| `OnRecordingStarted` | Fires when recording begins |
| `OnRecordingFinished(CarlRecording)` | Fires when recording stops |

**Methods:**

```csharp
void StartRecording(ulong maxSeconds = ulong.MaxValue);
CarlRecording StopRecording();
CarlExample CreateExample(CarlRecording recording, double startTimestamp, double endTimestamp);
```

**Properties:**

| Property | Type | Description |
|----------|------|-------------|
| `IsRecording` | `bool` | Whether recording is in progress |
| `ElapsedTime` | `double` | Seconds since `StartRecording()` |

## Programmatic Usage

Beyond the component workflow, you can use CARL's managed wrapper classes directly in C# scripts.

### Building Definitions in Code

```csharp
using Carl;

// Create a session (usually from CarlSessionManager.Instance.Session)
CarlSession session = CarlSessionManager.Instance.Session;

// Record input
CarlInProgressRecording inProgress = new CarlInProgressRecording(maxSeconds: 10);
// ... add samples each frame via inProgress.AddSample(ref sample) ...
CarlRecording recording = inProgress.Finish();

// Create examples with known action boundaries
CarlExample example1 = CarlExample.Create(recording, startTimestamp: 1.0, endTimestamp: 2.5);
CarlExample example2 = CarlExample.Create(recording, startTimestamp: 4.0, endTimestamp: 5.3);

// Build a definition
CarlDefinition definition = new CarlDefinition(ActionType.RightHandGesture);
definition.AddExample(recording, 1.0, 2.5);
definition.AddExample(recording, 4.0, 5.3);

// Save for later use
definition.SaveToFile(Application.persistentDataPath + "/my_gesture.carl");

// Clean up
example1.Dispose();
example2.Dispose();
recording.Dispose();
```

### Creating Recognizers in Code

```csharp
// Load a definition
CarlDefinition definition = CarlDefinition.LoadFromFile(path);

// Create a recognizer asynchronously
session.CreateRecognizerAsync(definition, recognizer => {
    // This runs on the main thread (after TickCallbacks)
    double score = recognizer.CurrentScore; // thread-safe
    recognizer.SetSensitivity(1.2);

    // Get canonical recording for visualization
    CarlRecordingInspector inspector = recognizer.GetCanonicalRecordingInspector();
    // ... use inspector ...
    inspector.Dispose();
});
```

## Action Types in Unity

The C# `ActionType` enum maps directly to the C++ `carl::action::ActionType`:

```csharp
public enum ActionType : int {
    LeftHandPose = 0,
    LeftHandGesture = 1,
    RightHandPose = 2,
    RightHandGesture = 3,
    TwoHandGesture = 4,
    LeftControllerGesture = 5,
    RightControllerGesture = 6,
    TwoControllerGesture = 7,
    LeftWristTrajectory = 8,
    RightWristTrajectory = 9,
    LeftHandShape = 10,
    RightHandShape = 11,
    Custom = -1,
}
```

See [Action Types](../../Documentation/ActionTypes.md) for detailed descriptions and selection guidance.

## Sensitivity Tuning

### On a Definition

```csharp
CarlDefinition definition = new CarlDefinition(ActionType.RightHandGesture);
definition.DefaultSensitivity = 0.8; // harder to trigger
```

### On a Recognizer

```csharp
recognizer.SetSensitivity(1.2); // easier to trigger (overrides definition default)
```

### On CarlGestureRecognizer (Component)

Set the `sensitivity` field in the Inspector or via script. Set to `-1` to use the definition's default.

Adjust `activationThreshold` and `deactivationThreshold` to control when the Unity Events fire.

## Auto-Trimming

Create precisely-trimmed examples using an existing recognizer:

```csharp
// Record a new performance
CarlRecording newRecording = recorder.StopRecording();

// Auto-trim using the recognizer's knowledge
CarlExample trimmed = CarlExample.CreateAutoTrimmed(recognizer, newRecording);

// Use the trimmed example
definition.AddExample(trimmed.GetRecording(), trimmed.StartTimestamp, trimmed.EndTimestamp);

trimmed.Dispose();
newRecording.Dispose();
```

This is useful for iteratively improving definitions — start with rough examples, then use auto-trimming to refine.

## Counterexamples

Add negative examples to help CARL distinguish similar gestures:

```csharp
// Record the motion that's being falsely recognized
CarlRecording confusingRecording = recorder.StopRecording();

// Add as a counterexample
definition.AddCounterexample(confusingRecording, startTimestamp, endTimestamp);

confusingRecording.Dispose();
```

## CarlDefinitionAsset

`CarlDefinitionAsset` is a `ScriptableObject` that stores a serialized gesture definition in the Unity asset system.

**Creating:**
- Right-click in Project → **Create → CARL → Definition Asset**

**Inspector:**
- Shows the `ActionType`, description, and data size
- **Import from File** button to load a `.carl` file
- **Export to File** button to save the current data

**Programmatic access:**

```csharp
// Load the definition from the asset
CarlDefinition definition = definitionAsset.Load();
// ... use it ...
definition.Dispose();

// Import from file
definitionAsset.ImportFromFile("/path/to/gesture.carl");

// Store a definition into the asset
definitionAsset.SetFromDefinition(definition, ActionType.RightHandGesture, "Wave Gesture");
```

## API Quick Reference

### Unity Components (MonoBehaviours)

| Class | Description |
|-------|-------------|
| `CarlSessionManager` | Singleton; owns the `CarlSession`; calls `TickCallbacks()` each frame |
| `CarlXRInputProvider` | Reads XR hardware; feeds InputSamples to the session |
| `CarlGestureRecognizer` | High-level gesture detection with hysteresis and UnityEvents |
| `CarlGestureRecorder` | Records gestures at runtime |

### Managed Wrapper Classes (IDisposable)

| Class | Description |
|-------|-------------|
| `CarlSession` | Runtime context; receives input; creates recognizers |
| `CarlInProgressRecording` | Mutable rolling buffer of InputSamples |
| `CarlRecording` | Immutable finalized recording |
| `CarlRecordingInspector` | Timestamp-based access to recording samples |
| `CarlExample` | Recording + start/end timestamps |
| `CarlDefinition` | ActionType + examples + counterexamples + sensitivity |
| `CarlRecognizer` | Scores live input against a definition |

### ScriptableObject

| Class | Description |
|-------|-------------|
| `CarlDefinitionAsset` | Unity asset wrapper for a serialized definition |

### Utility

| Class | Description |
|-------|-------------|
| `CarlCoordinateConversion` | Unity ↔ CARL coordinate system conversion (static) |

**Important:** All managed wrapper classes implement `IDisposable`. Always call `Dispose()` or use `using` blocks to avoid native memory leaks.

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
ninja carl_shared_library
# Copy libcarl.so to Targets/unity/Runtime/Plugins/Android/arm64-v8a/
```

## License

MIT — see [LICENSE.md](LICENSE.md).
