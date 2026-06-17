# FAQ

Common questions about using CARL.

## General

### What coordinate system does CARL use?

CARL uses a **right-handed Y-up** coordinate system, consistent with the [OpenXR](https://registry.khronos.org/OpenXR/specs/1.0/html/xrspec.html) convention. All 3D poses in `InputSample` must be in this coordinate system.

If you're working in Unity (which uses left-handed Y-up), the `CarlXRInputProvider` handles conversion automatically. For manual conversion: negate the Z component of positions, and negate the Z and W components of quaternions.

### How does CARL work internally?

CARL uses Dynamic Time Warping (DTW) to compare live input against stored examples. Each `ActionType` maps to a set of descriptors that extract relevant features from `InputSample`s. DTW handles variations in speed and timing, while descriptors handle variations in position and orientation by working in head-relative (egocentric) coordinate frames.

You don't need to understand DTW or descriptors to use CARL effectively — the library handles all of this automatically based on your choice of `ActionType`.

### Is CARL machine learning?

No. CARL uses classical (non-ML) techniques. There are no neural networks, no training phases, and no large datasets. A gesture can be defined from as few as 1–3 examples. This makes CARL suitable for on-device, real-time gesture definition by end users.

## Input

### What framerate does CARL need?

CARL has no minimum framerate requirement. Samples can arrive at any frequency. However, higher framerates generally produce smoother recognition scores. 30+ fps is recommended; 60+ fps is ideal for responsive recognition.

### Do timestamps need to be absolute (wall clock)?

No. Timestamps only need to be consistent within a session — they must use the same time reference and be monotonically increasing. You can use any epoch as long as all timestamps share it.

### Which InputSample fields should I populate?

Only populate the fields relevant to your `ActionType`. See [Action Types](ActionTypes.md) for which fields each type uses. The `Timestamp` field is always required. Unpopulated optional fields are ignored.

## Definitions

### How many examples should I provide?

**3–5 examples** is a good starting point for most gestures. More examples help CARL understand the natural variation in how a gesture is performed. Fewer examples (even 1) can work for simple, distinctive gestures.

Include natural variation in your examples — perform the gesture at slightly different speeds, from slightly different positions, and with natural hand differences.

### Can I add examples to an existing definition?

Yes. Call `definition.addExample()` to add more examples at any time. You'll need to create a new `Recognizer` from the updated definition for the changes to take effect.

### What makes a good counterexample?

A good counterexample is a motion that is similar to the target gesture but semantically different. For example, if your gesture is a "thumbs up" and it gets confused with a "pointing" pose, record the pointing pose and add it as a counterexample.

Don't add unrelated motions (like waving) as counterexamples for a thumbs up — they won't help because CARL wouldn't confuse them in the first place.

## Recognition

### How should I choose the right ActionType?

See the [decision tree in Action Types](ActionTypes.md#choosing-the-right-action-type). Key questions:
1. What input device? (hand tracking vs. controller)
2. Static pose or dynamic gesture?
3. One hand or two?
4. Does wrist orientation matter? (HandPose vs. HandShape)

### How many recognizers can I run simultaneously?

There is no hard limit. Each recognizer adds processing work, but CARL is designed to be lightweight. On mobile XR hardware (Quest), 10–20 recognizers running simultaneously is typical. Performance depends on the complexity of the definitions (number of examples, action type).

### What does the score mean?

`currentScore()` returns a confidence value:
- **≥ 1.0** — strong match; the action is being recognized
- **0.5–1.0** — partial match; the action may be in progress
- **≈ 0.0** — no match
- Can go slightly negative or above 1.0

For simple use, treat ≥ 1.0 as "recognized" and < 1.0 as "not recognized." The `whenRecognitionChangedSignal` handles this thresholding for you.

### Can I adjust recognition thresholds?

Yes. Use `Definition::DefaultSensitivity` or `Recognizer::setSensitivity()`. Higher values (> 1.0) make recognition easier to trigger; lower values (< 1.0) require closer matches.

In Unity, the `CarlGestureRecognizer` component provides `activationThreshold` and `deactivationThreshold` for hysteresis-based detection.

## Threading

### Is CARL thread-safe?

`Session` processes input on a background thread by default. Key threading rules:
- Call `addInput()` from your main thread
- Call `tickCallbacks()` from your main thread to dispatch signal handlers
- `currentScore()` is thread-safe and can be polled from any thread
- Signal handlers fire on the thread that calls `tickCallbacks()`

### When should I use single-threaded mode?

Use `Session{true}` (single-threaded) for:
- Debugging — all processing happens inline, making it easier to trace
- Environments without threading support (e.g., some WebAssembly contexts)
- Unit tests

For production XR applications, use the default multi-threaded mode.

## Migration

### I see "DescriptorType" in old code. What happened?

`DescriptorType` has been renamed to `ActionType`. The values and semantics are identical — only the name changed. Update your code to use `ActionType` and `getDescriptorType()` (which returns an `ActionType` value).

### Old docs say "left-handed Y-up." Is that still correct?

No. CARL uses **right-handed Y-up** (OpenXR convention). Earlier documentation incorrectly stated left-handed Y-up. If your code was written based on the old docs and works correctly, your data was likely already in the correct right-handed coordinate system. The only change is in the documentation, not the library behavior.
