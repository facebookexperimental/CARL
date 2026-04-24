# Tutorial: Defining Actions

This tutorial walks through the process of recording input, creating examples, and building a definition that CARL can use for recognition.

**Prerequisites:** Familiarity with [Concepts](../Concepts.md) and [Action Types](../ActionTypes.md).

## Step 1: Setting Up InputSamples

An `InputSample` is a snapshot of XR input state at a single instant. Which fields you populate depends on the `ActionType` you plan to use.

### For Hand Tracking Actions

For `LeftHandGesture`, `RightHandPose`, etc., populate hand joint poses:

```cpp
#include <carl/Carl.h>

using namespace carl;

InputSample createHandSample(double timestamp,
                              const TransformT& hmdPose,
                              const std::array<TransformT, 26>& rightJoints,
                              const TransformT& rightWrist) {
    InputSample sample{};
    sample.Timestamp = timestamp;
    sample.HmdPose = hmdPose;
    sample.RightWristPose = rightWrist;
    sample.RightHandJointPoses = rightJoints;
    return sample;
}
```

### For Controller Actions

For `LeftControllerGesture`, `RightControllerGesture`, etc., populate controller input:

```cpp
InputSample createControllerSample(double timestamp,
                                    const TransformT& hmdPose,
                                    const TransformT& rightWrist,
                                    float trigger, float squeeze,
                                    float thumbX, float thumbY) {
    InputSample sample{};
    sample.Timestamp = timestamp;
    sample.HmdPose = hmdPose;
    sample.RightWristPose = rightWrist;
    sample.RightControllerInput = std::array<NumberT, 7>{
        0.f,      // primary click
        0.f,      // secondary click
        thumbX,   // thumbstick X
        thumbY,   // thumbstick Y
        0.f,      // thumbstick click
        squeeze,  // squeeze value
        trigger   // trigger value
    };
    return sample;
}
```

### What to Populate

| ActionType | Required Fields |
|------------|----------------|
| Left/RightHandPose | Hand joint poses, wrist pose, HMD pose |
| Left/RightHandGesture | Hand joint poses, wrist pose, HMD pose |
| Left/RightHandShape | Hand joint poses only |
| TwoHandGesture | Both hands' joint + wrist poses, HMD pose |
| Left/RightControllerGesture | Controller input, wrist pose, HMD pose |
| TwoControllerGesture | Both controllers' input + wrist poses, HMD pose |
| Left/RightWristTrajectory | Wrist pose, HMD pose |

## Step 2: Recording Input

Use an `InProgressRecording` to accumulate InputSamples as they arrive:

```cpp
using namespace carl::action;

// Create a recording with a 10-second rolling window
InProgressRecording inProgress{/* maxSeconds */ 10};

// In your input loop (e.g., 60 fps):
void onNewFrame(double timestamp, /* ... input data ... */) {
    InputSample sample = createHandSample(timestamp, /* ... */);
    inProgress.addSample(std::move(sample));
}
```

The `maxSeconds` parameter limits how much history is kept. For gesture recording, 10 seconds is typically sufficient. Use the default constructor (`InProgressRecording{}`) for unlimited recording.

When you are done recording, finalize into an immutable `Recording`:

```cpp
Recording recording{std::move(inProgress)};
// inProgress is now consumed and cannot be used
```

## Step 3: Creating an Example

An `Example` marks when an action occurred within a `Recording` by specifying start and end timestamps:

```cpp
// The action occurred between 2.5 and 3.8 seconds
Example example{recording, /* startTimestamp */ 2.5, /* endTimestamp */ 3.8};
```

### Choosing Timestamps

You need to determine which portion of the recording contains the action. There are several approaches:

1. **Manual specification** — if you know when the user started and stopped the action (e.g., from a UI button press), use those timestamps directly.

2. **Using RecordingInspector** — scrub through the recording to find the right boundaries:

   ```cpp
   RecordingInspector inspector = recording.getInspector();
   double start = inspector.startTimestamp(); // first sample
   double end = inspector.endTimestamp();     // last sample

   // Inspect the input at any timestamp to verify boundaries
   const InputSample& sample = inspector.inspect(2.5);
   ```

3. **Auto-trimming** — if you already have a working recognizer, use it to automatically find the best boundaries. See [Tuning and Debugging](TuningAndDebugging.md#auto-trimming).

## Step 4: Building a Definition

A `Definition` combines an `ActionType` with one or more examples:

```cpp
Definition definition{ActionType::RightHandGesture};

// Add examples (3-5 recommended for good recognition)
definition.addExample(std::move(example1));
definition.addExample(std::move(example2));
definition.addExample(std::move(example3));
```

### Adding Counterexamples

If similar actions might be confused, add counterexamples to help CARL distinguish them:

```cpp
// Record a similar-but-different action
definition.addCounterexample(std::move(similarButDifferent));
```

See [Tuning and Debugging](TuningAndDebugging.md#counterexamples) for guidance.

### Setting Sensitivity

The default sensitivity is 1.0. Adjust it to make the action easier or harder to trigger:

```cpp
definition.DefaultSensitivity = 0.8; // harder to trigger (requires closer match)
definition.DefaultSensitivity = 1.2; // easier to trigger (accepts looser matches)
```

## Step 5: Saving the Definition

Save the definition for later use:

```cpp
#include <carl/utilities/FileSerialization.h>

carl::utilities::SerializeToFile(definition, "my_gesture.carl");
```

Load it back:

```cpp
auto loaded = carl::utilities::TryDeserializeFromFile<Definition>("my_gesture.carl");
if (loaded.has_value()) {
    Definition& definition = loaded.value();
    // Use the definition...
}
```

See [Serialization](../Serialization.md) for more details on save/load options.

## Complete Example

Putting it all together — recording a right-hand gesture and saving the definition:

```cpp
#include <carl/Carl.h>
#include <carl/utilities/FileSerialization.h>

using namespace carl;
using namespace carl::action;

void recordAndSave() {
    // 1. Record input
    InProgressRecording inProgress{10};
    // ... add samples in your input loop ...
    Recording recording{std::move(inProgress)};

    // 2. Create examples with known boundaries
    Example example1{recording, 1.0, 2.5};
    Example example2{recording, 4.0, 5.3};
    Example example3{recording, 7.0, 8.1};

    // 3. Build the definition
    Definition definition{ActionType::RightHandGesture};
    definition.addExample(std::move(example1));
    definition.addExample(std::move(example2));
    definition.addExample(std::move(example3));

    // 4. Save it
    carl::utilities::SerializeToFile(definition, "my_gesture.carl");
}
```

## Next Steps

- [Recognizing Actions](RecognizingActions.md) — use this definition for live recognition
- [Tuning and Debugging](TuningAndDebugging.md) — improve recognition accuracy
- [Serialization](../Serialization.md) — more on saving and loading
