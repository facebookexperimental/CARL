# CARL Unity Integration — Documentation

## Overview

The CARL Unity package provides gesture recognition for XR (VR/AR) applications using the **Classical Action Recognition Library**. CARL uses **Dynamic Time Warping (DTW)** to match gestures — no machine learning or large datasets required. Gestures are defined simply by performing them.

## Architecture

The package has three layers:

### 1. Native P/Invoke Layer (`Carl.Native`)
Raw bindings to the CARL C API (`carl.h`). You typically don't interact with this layer directly.

- `CarlNative` — `DllImport` declarations for all 42 C API functions.
- `CarlInteropTypes` — Marshaled structs matching the native memory layout (`CarlInputSampleInterop`, `CarlOptionalTransform`, etc.).
- `CarlBytes` — Helper for the two-call byte buffer pattern used by serialization functions.

### 2. Managed Wrapper Layer (`Carl`)
Idiomatic C# classes wrapping native handles with `IDisposable` lifecycle management.

| Class | Purpose |
|-------|---------|
| `CarlSession` | Central session managing input processing and recognizers |
| `CarlDefinition` | Defines an action with examples and counterexamples |
| `CarlRecognizer` | Evaluates live input against a definition, producing confidence scores |
| `CarlRecording` | A completed recording of InputSamples |
| `CarlInProgressRecording` | An in-progress recording accumulating samples |
| `CarlRecordingInspector` | Allows scrubbing through a recording by timestamp |
| `CarlExample` | An action occurrence within a recording (start/end timestamps) |

### 3. Unity Component Layer (`Carl`)
`MonoBehaviour` and `ScriptableObject` components for scene integration.

| Component | Purpose |
|-----------|---------|
| `CarlSessionManager` | Singleton that creates/manages the session and pumps callbacks |
| `CarlXRInputProvider` | Reads XR hand/controller/HMD input and feeds CARL each frame |
| `CarlGestureRecognizer` | Threshold-based gesture detection with UnityEvents |
| `CarlGestureRecorder` | Records InputSamples at runtime to create new gestures |
| `CarlDefinitionAsset` | ScriptableObject for storing gesture definitions as Unity assets |

## Quick Start

### 1. Set Up the Scene

Add these components to GameObjects in your scene:

```
[GameObject: CARL Manager]
  └─ CarlSessionManager (singleton, auto-creates session)
  └─ CarlXRInputProvider (reads XR input, feeds CARL)
```

### 2. Create a Definition Asset

1. In the Project window, right-click → **Create → CARL → Definition Asset**.
2. In the Inspector, click **Import from File** and select a `.carl` definition file.

### 3. Add a Gesture Recognizer

```
[GameObject: My Gesture]
  └─ CarlGestureRecognizer
       Definition Asset: [drag your asset here]
       Activation Threshold: 0.8
       Deactivation Threshold: 0.6
       OnGestureDetected: [wire up your handler]
       OnGestureLost: [wire up your handler]
```

### 4. Respond to Gestures

```csharp
public class MyGestureHandler : MonoBehaviour
{
    public void OnGestureDetected()
    {
        Debug.Log("Gesture detected!");
    }

    public void OnGestureLost()
    {
        Debug.Log("Gesture ended.");
    }
}
```

## Coordinate System

| System | Convention |
|--------|-----------|
| Unity | Left-handed, Y-up |
| CARL / OpenXR | Right-handed, Y-up |

**Conversion (handled automatically by `CarlXRInputProvider`):**
- Positions: negate Z component
- Quaternions: negate X and Y components

If you feed InputSamples manually, use `CarlCoordinateConversion`:

```csharp
var carlTransform = CarlCoordinateConversion.ToCarlTransform(unityPose);
var unityPose = CarlCoordinateConversion.ToUnityPose(carlTransform);
```

## Threading Model

- `CarlSession` processes gesture matching on a **background thread** by default.
- `CarlSessionManager.Update()` calls `TickCallbacks()` each frame, dispatching results (recognizer creation callbacks, etc.) to the **main thread**.
- `CarlRecognizer.CurrentScore` is **thread-safe** and can be polled at any time.
- Pass `singleThreaded: true` to `CarlSession` constructor (or check the box on `CarlSessionManager`) to process everything inline — useful for debugging.

## Action Types

| ActionType | Description |
|------------|-------------|
| `LeftHandPose` / `RightHandPose` | Static hand shape + wrist orientation |
| `LeftHandGesture` / `RightHandGesture` | Dynamic hand gesture (shape + motion) |
| `TwoHandGesture` | Gesture involving both hands |
| `LeftControllerGesture` / `RightControllerGesture` | VR controller gesture |
| `TwoControllerGesture` | Both controllers |
| `LeftWristTrajectory` / `RightWristTrajectory` | Wrist motion path only |
| `LeftHandShape` / `RightHandShape` | Hand shape without orientation |

## Recording Gestures at Runtime

Use `CarlGestureRecorder` to capture gestures:

```csharp
var recorder = GetComponent<CarlGestureRecorder>();

// Start recording
recorder.StartRecording();

// ... user performs gesture ...

// Stop recording
CarlRecording recording = recorder.StopRecording();

// Create an example from the full recording
using var example = CarlExample.Create(recording,
    recording.GetInspector().StartTimestamp,
    recording.GetInspector().EndTimestamp);

// Save for later use
example.SaveToFile("my_gesture.carl_example");
```

## API Reference

### CarlSession

```csharp
var session = new CarlSession(singleThreaded: false);
session.SetLogger(msg => Debug.Log(msg));
session.AddInputSample(ref sample);
session.TickCallbacks(); // call every frame
session.CreateRecognizerAsync(definition, recognizer => { ... });
session.Dispose();
```

### CarlDefinition

```csharp
var def = new CarlDefinition(ActionType.RightHandGesture);
def.AddExample(recording, startTime, endTime);
def.DefaultSensitivity = 5.0;
byte[] bytes = def.Serialize();
def.SaveToFile("gesture.carl");

var loaded = CarlDefinition.LoadFromFile("gesture.carl");
var deserialized = CarlDefinition.Deserialize(bytes);
```

### CarlRecognizer

```csharp
session.CreateRecognizerAsync(definition, recognizer =>
{
    recognizer.SetSensitivity(3.0);
    double score = recognizer.CurrentScore; // 0.0 - 1.0+
    recognizer.Dispose(); // async deletion
});
```

## Building Native Plugins

### Windows (x64)
```bash
mkdir build && cd build
cmake ../CARL
cmake --build . --config Release
```
Copy `carl.dll` to `Runtime/Plugins/x86_64/`.

### Android / Quest (arm64-v8a)
```bash
mkdir build-android && cd build-android
cmake -G Ninja ../CARL \
  -DCMAKE_SYSTEM_NAME=Android \
  -DCMAKE_ANDROID_ARCH_ABI="arm64-v8a" \
  -DCMAKE_ANDROID_NDK=/path/to/ndk \
  -DCMAKE_ANDROID_STL_TYPE=c++_static \
  -DCMAKE_BUILD_TYPE=Release
ninja carl
```
Copy `libcarl.so` to `Runtime/Plugins/Android/arm64-v8a/`.

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `DllNotFoundException: carl` | Native plugin binary not found. See build instructions above. |
| Struct size mismatch error on startup | ABI incompatibility — rebuild native plugin with matching compiler/settings. |
| `CarlGestureRecognizer` never fires | Check that `CarlXRInputProvider` is in the scene and XR hand tracking is active. |
| Score is always 0.0 | Verify the definition has examples and the correct `ActionType` for your input. |
