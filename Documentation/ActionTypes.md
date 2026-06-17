# Action Types

CARL supports 12 built-in action types and a custom type. Each action type determines which fields from an `InputSample` are used for recognition, how those fields are compared between frames, and whether the action is static (a pose) or dynamic (a gesture over time).

This page describes each action type and provides guidance on choosing the right one.

## Hand Tracking Types

These action types use hand joint data from `InputSample.LeftHandJointPoses` / `RightHandJointPoses` and wrist poses.

### LeftHandPose / RightHandPose

**Static single-hand pose.** Recognizes a hand shape combined with wrist orientation — for example, a thumbs-up, a pointing gesture, or an "OK" sign at a specific angle.

- **Input used:** 26 hand joint poses + wrist orientation (relative to HMD)
- **Static:** Recognized as long as the pose is held; no motion component
- **Use when:** The gesture is defined by both finger configuration *and* wrist angle

### LeftHandGesture / RightHandGesture

**Dynamic single-hand gesture.** Recognizes a hand motion over time — for example, a wave, a swipe, or a snap.

- **Input used:** 26 hand joint poses + wrist pose + wrist motion + wrist rotation (all relative to HMD)
- **Dynamic:** Requires the user to perform the motion; not just hold a pose
- **Use when:** The gesture involves motion through space and/or changing hand shape over time

### TwoHandGesture

**Dynamic two-handed gesture.** Recognizes coordinated motion of both hands — for example, a clap, a pulling-apart motion, or sign language gestures.

- **Input used:** All data from both hands + relative position between hands
- **Dynamic:** Both hands' motion over time, including their spatial relationship
- **Use when:** The gesture requires coordination between both hands

### LeftHandShape / RightHandShape

**Static finger configuration only.** Recognizes only the shape of the fingers, ignoring wrist orientation entirely — for example, an open palm vs. a fist, regardless of which direction the hand is pointing.

- **Input used:** 26 hand joint poses only (wrist-relative, no orientation)
- **Static:** Recognized as long as the shape is held
- **Use when:** You care about finger configuration but not hand orientation (e.g., distinguishing open vs. closed hand regardless of where the user is pointing)

### How HandShape differs from HandPose

`HandPose` includes wrist orientation, so a thumbs-up pointing left is different from a thumbs-up pointing right. `HandShape` ignores orientation, so both would match. Choose `HandShape` when orientation does not matter.

## Controller Types

These action types use controller button/axis data from `InputSample.LeftControllerInput` / `RightControllerInput` along with wrist motion.

### Controller Input Channels

Each controller provides 7 input channels:

| Index | Channel | Range |
|-------|---------|-------|
| 0 | Primary click | 0.0 or 1.0 |
| 1 | Secondary click | 0.0 or 1.0 |
| 2 | Thumbstick X | -1.0 to 1.0 |
| 3 | Thumbstick Y | -1.0 to 1.0 |
| 4 | Thumbstick click | 0.0 or 1.0 |
| 5 | Squeeze value | 0.0 to 1.0 |
| 6 | Trigger value | 0.0 to 1.0 |

### LeftControllerGesture / RightControllerGesture

**Dynamic single-controller gesture.** Recognizes a combination of button/axis input and controller motion — for example, pulling the trigger while sweeping right, or a thumbstick circle.

- **Input used:** 7 controller channels + wrist trajectory + wrist orientation + wrist rotation
- **Dynamic:** Motion over time with controller state
- **Use when:** The gesture involves controller buttons/axes, optionally combined with hand motion

### TwoControllerGesture

**Dynamic two-controller gesture.** Recognizes coordinated motion of both controllers.

- **Input used:** All data from both controllers + relative position
- **Dynamic:** Both controllers' motion and input over time
- **Use when:** The gesture requires coordination between both controllers

## Wrist Trajectory Types

### LeftWristTrajectory / RightWristTrajectory

**Dynamic wrist path.** Recognizes the path a wrist takes through space, ignoring finger configuration entirely — for example, a circular motion, a figure-eight, or a throw.

- **Input used:** Wrist position over time + timing information (relative to HMD)
- **Dynamic:** Path through space with temporal matching
- **Use when:** Only the motion path matters, not what the fingers are doing

## Custom Type

### Custom

**User-defined action type.** Define your own descriptor logic by implementing the required operations and registering with `Session::enableCustomActionType<T>()`.

You must provide:
- `TryCreate` — extract a descriptor from current and prior InputSamples
- `Distance` — compute distance between two descriptors
- `Lerp` — interpolate between two descriptors
- `CalculateTuning` — compute per-dimension tuning weights from training examples

See the [C++ API Reference](CppApiReference.md#session) for the full `TemplatedCustomActionTypeOperations` interface.

## Choosing the Right Action Type

Use this decision tree to select an action type:

```
What input device?
├── Hand tracking
│   ├── Is it a static pose (no motion)?
│   │   ├── Does wrist orientation matter?
│   │   │   ├── Yes → LeftHandPose / RightHandPose
│   │   │   └── No  → LeftHandShape / RightHandShape
│   │   └──
│   ├── Is it a dynamic gesture (motion over time)?
│   │   ├── One hand  → LeftHandGesture / RightHandGesture
│   │   └── Two hands → TwoHandGesture
│   └── Only wrist motion matters (ignore fingers)?
│       └── LeftWristTrajectory / RightWristTrajectory
│
├── Controller
│   ├── One controller  → LeftControllerGesture / RightControllerGesture
│   └── Two controllers → TwoControllerGesture
│
└── Something else → Custom
```

### Tips

- **Start with the most specific type.** If wrist orientation matters, use `HandPose` not `HandShape`. You can always switch to a broader type if recognition is too strict.
- **Poses vs. gestures:** Poses are recognized continuously while held. Gestures fire once when the motion is completed.
- **Counterexamples help:** If two similar gestures are being confused, add counterexamples rather than switching action types. See [Tuning and Debugging](Tutorials/TuningAndDebugging.md).
- **Multiple examples:** Provide 3–5 examples of each gesture for best results, with natural variation in speed and style.
