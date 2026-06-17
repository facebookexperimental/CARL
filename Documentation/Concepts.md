# CARL Concepts

This document describes the core concepts of CARL and how they fit together. After reading this, you will understand the data model that drives action recognition and be ready to follow the [tutorials](Tutorials/).

## Introduction

CARL is inspired by the [$-family of recognizers](https://depts.washington.edu/acelab/proj/dollar/impact.html) — a family of algorithms that tackle action recognition using classical techniques rather than machine learning. CARL's guiding principle is that **actions are best defined by simply doing them**. Perform a gesture a few times, and CARL learns to recognize it. No training data, no statistical models, no specialized expertise.

This makes it possible for any user who can perform a motion to create action definitions on the fly, opening the door to personalized gesture recognition in XR applications.

## Coordinate System

CARL uses a **right-handed Y-up** coordinate system, consistent with the [OpenXR](https://registry.khronos.org/OpenXR/specs/1.0/html/xrspec.html) convention:

- **+X** points right
- **+Y** points up
- **+Z** points toward the user (out of the screen)

All 3D poses provided to CARL must be expressed in this coordinate system. If you are working in Unity (which uses left-handed Y-up), the `CarlXRInputProvider` handles conversion automatically. For manual conversion, negate the Z component of positions and the Z and W components of quaternions. See [Unity Integration](../Targets/unity/README.md) for details.

## InputSample

An `InputSample` is CARL's atomic unit of input: a snapshot of all relevant XR state at a single instant. It contains:

| Field | Type | Description |
|-------|------|-------------|
| `Timestamp` | `double` | Time in seconds (required) |
| `HmdPose` | optional | Head-mounted display 6DOF pose (position + orientation) |
| `LeftWristPose` | optional | Left wrist 6DOF pose |
| `RightWristPose` | optional | Right wrist 6DOF pose |
| `LeftHandJointPoses` | optional | 26 left-hand joint poses (see below) |
| `RightHandJointPoses` | optional | 26 right-hand joint poses |
| `LeftControllerInput` | optional | 7 left controller input channels |
| `RightControllerInput` | optional | 7 right controller input channels |

**Timestamp** is the only required field. Everything else is optional because not all data may be available at all times (e.g., a hand may not be tracked). Which fields you populate depends on the [ActionType](ActionTypes.md) you intend to use.

**Hand joints** follow the [OpenXR XR_EXT_hand_tracking](https://registry.khronos.org/OpenXR/specs/1.0/html/xrspec.html#XR_EXT_hand_tracking) layout: 26 joints per hand, from Palm (0) through Little Fingertip (25).

**Controller input channels** map to standard XR controller inputs:

| Index | Channel |
|-------|---------|
| 0 | Primary click |
| 1 | Secondary click |
| 2 | Thumbstick X |
| 3 | Thumbstick Y |
| 4 | Thumbstick click |
| 5 | Squeeze value |
| 6 | Trigger value |

## Recording

A `Recording` is a continuous, time-ordered sequence of `InputSample`s. Recordings have no fixed framerate requirement — samples can arrive at any frequency. The only requirement is that all timestamps within a recording use the same time reference.

You build a `Recording` through an `InProgressRecording`, which acts as a rolling buffer:

```cpp
InProgressRecording inProgress{/* maxSeconds */ 10}; // keeps last 10 seconds
inProgress.addSample(sample1);
inProgress.addSample(sample2);
// ...
Recording recording{std::move(inProgress)}; // finalize
```

The `maxSeconds` parameter limits how much history is retained. When samples older than `maxSeconds` are pushed out, they are discarded. Use the default constructor for an unlimited buffer.

### RecordingInspector

A `RecordingInspector` provides efficient, timestamp-based access to a recording's samples. Given a timestamp, it returns the closest `InputSample`:

```cpp
RecordingInspector inspector = recording.getInspector();
double start = inspector.startTimestamp();
double end = inspector.endTimestamp();
const InputSample& sample = inspector.inspect(someTimestamp);
```

A `RecordingInspector` is a non-owning view — the underlying `Recording` must outlive it.

## Example

An `Example` pairs a `Recording` with a start and end timestamp that mark when an action occurred within that recording:

```cpp
Example example{recording, startTimestamp, endTimestamp};
```

The recording typically contains more data than just the action itself (e.g., idle time before and after). The timestamps identify the relevant portion.

An `Example` contains no information about *what* action occurred — it only says that *something* happened between those timestamps. The semantic meaning comes from the `Definition` it is added to.

## Counterexample

A counterexample has the same structure as an `Example` — a recording with start/end timestamps — but serves the opposite purpose. While examples teach CARL what an action looks like, counterexamples teach it what the action does **not** look like. This is useful for disambiguating similar gestures (e.g., distinguishing a wave from a swipe).

Add counterexamples to a `Definition` via `addCounterexample()`.

## Definition

A `Definition` brings together:

1. An **ActionType** — which kind of action this defines (see [Action Types](ActionTypes.md))
2. One or more **examples** — recordings of the action being performed
3. Zero or more **counterexamples** — recordings of similar-but-different actions
4. A **DefaultSensitivity** — how easily the action triggers (default: 1.0)

```cpp
Definition definition{ActionType::RightHandGesture};
definition.addExample(example1);
definition.addExample(example2);
definition.addCounterexample(similarButWrong);
definition.DefaultSensitivity = 0.8;
```

Definitions are serializable and can be saved to and loaded from `.carl` files. See [Serialization](Serialization.md).

## ActionType

CARL supports 12 built-in action types plus custom types. Each action type determines which `InputSample` fields are used for recognition and how they are compared:

| ActionType | Input Used |
|------------|------------|
| `LeftHandPose` / `RightHandPose` | Hand shape + wrist orientation (static) |
| `LeftHandGesture` / `RightHandGesture` | Full hand + wrist motion (dynamic) |
| `TwoHandGesture` | Both hands + relative position |
| `LeftHandShape` / `RightHandShape` | Finger configuration only (static) |
| `LeftControllerGesture` / `RightControllerGesture` | Controller buttons/axes + wrist motion |
| `TwoControllerGesture` | Both controllers |
| `LeftWristTrajectory` / `RightWristTrajectory` | Wrist path through space (ignoring fingers) |
| `Custom` | User-defined via `Session::enableCustomActionType()` |

For detailed descriptions and guidance on choosing the right type, see [Action Types](ActionTypes.md).

## Session

A `Session` is the runtime context for recognition. It receives `InputSample`s and hosts `Recognizer`s:

```cpp
Session session{/* singleThreaded */ false};

// Feed input each frame
session.addInput(liveSample);

// Dispatch callbacks on the calling thread
arcana::cancellation token{};
session.tickCallbacks(token);
```

By default, a `Session` processes recognition on a background thread. Pass `singleThreaded = true` to run everything on the calling thread (useful for debugging or single-threaded environments).

Use `setLogger()` to receive diagnostic messages:

```cpp
session.setLogger([](std::string msg) { std::cout << msg << std::endl; });
```

## Recognizer

A `Recognizer` watches a `Session` and scores incoming input against a `Definition`. It is the final piece that produces recognition results:

```cpp
Recognizer recognizer{session, definition};
```

### Confidence Score

`currentScore()` returns a `double` indicating how confident CARL is that the defined action is happening:

- **≥ 1.0** — strong match; the action is recognized
- **≈ 0.0** — no match
- Values can exceed 1.0 or go slightly below 0.0

For simple use cases, clamp to [0, 1] and treat as a confidence level.

### Recognition Events

Subscribe to `whenRecognitionChangedSignal` for discrete start/stop events:

```cpp
auto ticket = recognizer.whenRecognitionChangedSignal.addHandler(
    [](bool recognized) {
        if (recognized) { /* action started */ }
        else            { /* action ended   */ }
    });
```

The returned `ticket` controls the subscription lifetime — the handler is unsubscribed when the ticket is destroyed.

### Sensitivity

Each recognizer inherits its definition's `DefaultSensitivity` but can override it:

```cpp
recognizer.setSensitivity(1.2); // easier to trigger
```

Higher sensitivity means the action triggers more easily (lower score threshold).

### Additional Features

- **`createAutoTrimmedExample(recording)`** — uses the recognizer's knowledge to automatically determine optimal start/end timestamps for a new example from a raw recording
- **`analyzeRecording(recording, ostream)`** — outputs diagnostic information about how well a recording matches the definition
- **`getCanonicalRecordingInspector()`** — returns an inspector for the recognizer's internal canonical representation

## Lifecycle Overview

```
InputSample ─► InProgressRecording ─► Recording ─► Example ─► Definition
                                                          └─► Counterexample

Session.addInput(sample) ──────────────────────────► Recognizer(session, definition)
                                                          │
                                                          ├─ currentScore() → double
                                                          └─ whenRecognitionChangedSignal → bool
```

1. Collect `InputSample`s into an `InProgressRecording`
2. Finalize into a `Recording`
3. Mark action boundaries to create `Example`s
4. Bundle examples into a `Definition` with an `ActionType`
5. Create a `Session` and feed live `InputSample`s
6. Attach a `Recognizer` with your `Definition` and listen for recognition events
