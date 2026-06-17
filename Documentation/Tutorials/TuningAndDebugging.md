# Tutorial: Tuning and Debugging

This tutorial covers techniques for improving recognition quality: adjusting sensitivity, using counterexamples, auto-trimming examples, and diagnosing issues.

**Prerequisites:** Familiarity with [Defining Actions](DefiningActions.md) and [Recognizing Actions](RecognizingActions.md).

## Sensitivity

Sensitivity controls how easily an action triggers. It can be set at two levels:

### Definition-Level (Default)

```cpp
Definition definition{ActionType::RightHandGesture};
definition.DefaultSensitivity = 1.0; // default
```

All recognizers created from this definition inherit this value.

### Recognizer-Level (Override)

```cpp
Recognizer recognizer{session, definition};
recognizer.setSensitivity(1.2); // override for this instance
```

### How Sensitivity Works

| Value | Effect |
|-------|--------|
| < 1.0 | Harder to trigger — requires a closer match to the examples |
| 1.0 | Default behavior |
| > 1.0 | Easier to trigger — accepts looser matches |

**Tips:**
- Start with the default (1.0) and adjust based on testing
- If the action never triggers, try increasing sensitivity (e.g., 1.2–1.5)
- If the action triggers on unrelated motions, try decreasing sensitivity (e.g., 0.7–0.9)
- Small adjustments (0.1–0.2 increments) are usually sufficient

## Counterexamples

Counterexamples teach CARL what an action does **not** look like. They are especially useful when two similar actions are being confused.

### When to Use Counterexamples

- A wave gesture is being confused with a swipe
- A thumbs-up is triggering on a pointing gesture
- Any time the recognizer fires on motions that look similar but are semantically different

### How to Add Counterexamples

Record the confusing motion (the one that shouldn't trigger), create an example from it, and add it as a counterexample:

```cpp
// Record the motion that's being falsely recognized
Recording confusingRecording{std::move(inProgress)};
Example confusingExample{confusingRecording, startTime, endTime};

// Add as counterexample to the definition
definition.addCounterexample(std::move(confusingExample));
```

### Tips for Effective Counterexamples

- **Record naturally** — perform the confusing motion exactly as it happens in practice
- **Be specific** — add counterexamples for the specific motions causing false positives, not unrelated motions
- **Don't overdo it** — 1–3 counterexamples per confusing motion is usually enough
- **Test after each addition** — verify that the counterexample fixes the confusion without breaking valid recognition

## Auto-Trimming

When recording a gesture, it can be difficult to identify exactly when the action started and ended. Auto-trimming solves this by using an existing recognizer to find the optimal boundaries.

### Bootstrap Workflow

1. Create an initial definition with roughly-trimmed examples
2. Build a recognizer from that definition
3. Record new performances of the gesture
4. Use auto-trimming to create precisely-trimmed examples
5. Build a better definition from the auto-trimmed examples

```cpp
// 1. Start with a rough definition
Definition roughDef{ActionType::RightHandGesture};
roughDef.addExample(roughExample1);
roughDef.addExample(roughExample2);

// 2. Create a recognizer
Session session{true};
Recognizer recognizer{session, roughDef};

// 3. Record a new performance
Recording newRecording{std::move(inProgress)};

// 4. Auto-trim to find optimal boundaries
Example trimmed = recognizer.createAutoTrimmedExample(newRecording);

// 5. Build a refined definition
Definition refinedDef{ActionType::RightHandGesture};
refinedDef.addExample(std::move(trimmed));
// ... add more auto-trimmed examples ...
```

Auto-trimming uses the recognizer's knowledge of the action to determine where the most characteristic portion of the recording is. This generally produces better example boundaries than manual trimming.

## Analyzing Recordings

Use `analyzeRecording()` to get diagnostic information about how well a recording matches a definition:

```cpp
recognizer.analyzeRecording(newRecording, std::cout);
```

This outputs detailed information to the provided stream, including per-example match costs. Use this to understand:
- Why a particular recording isn't being recognized
- Which examples in the definition are contributing most to (or detracting from) the match
- Whether the recording quality is sufficient

## Common Issues and Fixes

### Action Never Triggers

**Symptoms:** `currentScore()` stays near 0; `whenRecognitionChangedSignal` never fires `true`.

**Possible causes and fixes:**
1. **Wrong ActionType** — verify that the action type matches your input. For example, `HandPose` won't work if you're only providing wrist data without hand joints.
2. **Missing input fields** — ensure your `InputSample` populates all fields required by the ActionType (see [Action Types](../ActionTypes.md)).
3. **Too few examples** — add more examples (3–5 recommended) with natural variation.
4. **Sensitivity too low** — increase `DefaultSensitivity` or per-recognizer sensitivity.
5. **Coordinate system mismatch** — all poses must be in right-handed Y-up (OpenXR). See [Concepts](../Concepts.md#coordinate-system).

### Action Triggers Too Easily

**Symptoms:** High scores on unrelated motions; frequent false positives.

**Possible causes and fixes:**
1. **Sensitivity too high** — decrease sensitivity.
2. **Too few examples** — more examples help CARL understand the range of valid performances. Counterintuitively, more examples can make recognition more specific, not less.
3. **No counterexamples** — add counterexamples for the motions causing false positives.
4. **Wrong ActionType** — consider a more specific action type. For example, use `HandPose` instead of `HandShape` if wrist orientation matters.

### Two Gestures Interfere With Each Other

**Symptoms:** Both recognizers fire at the same time, or the wrong one fires.

**Possible causes and fixes:**
1. **Add counterexamples** — for each gesture, add examples of the other gesture as counterexamples.
2. **Use different ActionTypes** — if the gestures use different input modalities, assign them different ActionTypes.
3. **Improve examples** — ensure each gesture's examples are distinct and clearly demonstrate the intended motion.

### Recognition Feels Laggy

**Symptoms:** Score changes slowly; events fire noticeably after the action is complete.

**Possible causes and fixes:**
1. **Input framerate too low** — feed InputSamples more frequently. 30+ fps is recommended; 60+ is ideal.
2. **Multi-threaded callback delay** — ensure `tickCallbacks()` is called every frame on the main thread.
3. **Example timestamps too wide** — tighten the start/end timestamps on your examples to match just the core action, not idle time before and after. Try auto-trimming.

## Next Steps

- [C++ API Reference](../CppApiReference.md) — full method documentation
- [Action Types](../ActionTypes.md) — revisit action type selection
- [Serialization](../Serialization.md) — save tuned definitions
