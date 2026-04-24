# Auto-Trimming & Adaptive Definition: Health Notes

Observations from code review, April 2026. These should be addressed before relying on auto-trimming in production scenarios.

---

## Bugs

### 1. Timestamp boundary fallthrough in `createAutoTrimmedExample`

**File:** `Core/CARL/CARL/source/Recognizer.cpp`, lines 235-248

When no input samples exist before `firstDescriptorT` or after `lastDescriptorT`, `startT` and `endT` remain at `numeric_limits<NumberT>::lowest()` / `::max()`. The returned `Example` has absurd timestamp bounds spanning the float range. This happens when the gesture occupies the very beginning or end of the recording.

**Fix:** Default `startT`/`endT` to the recording's actual first/last sample timestamps rather than float extremes.

### 2. Copy-paste bug in commented-out handedness filtering

**File:** `Core/CARL/CARL/source/Recognizer.cpp`, lines 446-448

Both lines in the `RightHanded` branch say `sampleCopy.LeftHandJointPoses.reset()`. The second should be `sampleCopy.LeftWristPose.reset()`. Dead code currently, but will bite if uncommented.

### 3. Recognition state machine is pulse-only

**File:** `Core/CARL/CARL/source/Recognizer.cpp`, lines 598-612

```cpp
if (!m_recognition && score > 0.7)       // → set true, fire
else if (m_recognition)                   // → set false, fire (ALWAYS, regardless of score)
```

Recognition resets on the very next frame after firing, regardless of whether score is still above threshold. This means `m_recognition` is never sustained — it's a one-frame pulse. For real-time observation of ongoing recognition (e.g., deciding when a user's attempt is "done"), you'd need to observe `currentScore()` directly rather than relying on the recognition state.

Consider whether this should instead be:
```cpp
else if (m_recognition && score <= 0.7)   // → only reset when score actually drops
```

---

## Concerning TODOs

### 4. `EgocentricTemporalSpace::getPose()` returns Identity on degenerate geometry

**File:** `Core/CARL/Utilities/include/carl/descriptor/EgocentricTemporalSpace.h`, line 27

When origin ≈ ego position (wrist near head), z-axis degenerates and the function returns `TransformT::Identity()`. Descriptors computed during that frame will be in world space rather than egocentric space, creating discontinuities in the sequence. Users whose wrists pass near their face mid-gesture could get bad auto-trimming.

### 5. `getWristPose()` implicit behavior change

**File:** `Core/CARL/Utilities/include/carl/descriptor/DescriptorUtils.h`, line 69

The author's own comment: *"This implicit behavior change based on the presence or absence of joint data is a bug farm."* With joints present, it computes a pose from joint geometry; without, it falls back to the raw wrist pose. The descriptor sequences are qualitatively different. If hand tracking drops joints intermittently, descriptor sequences become heterogeneous.

### 6. Delta descriptor interaction with resampling

**File:** `Core/CARL/Utilities/include/carl/descriptor/SequenceOperations.h`, lines 32-33

Open TODO: *"Figure out if delta descriptors (rotation, translation, etc.) will play correctly with this."* `extendSequence()` uses bisection/interpolation to maintain consistent descriptor spacing. Delta descriptors measure frame-to-frame change, and it's unclear whether interpolated intermediate samples produce correct delta values. Auto-trimming builds its descriptor sequence through `extendSequence()`, so this directly affects trim boundary accuracy.

### 7. Tuning resampling order

**File:** `Core/CARL/CARL/source/Recognizer.cpp`, lines 168-170

TODO: *"Figure out why the tuning resampling negatively impacts recognition, then substitute the following for the above."* Currently, templates are initialized before tuning is calculated — but ideally tuning should be computed first so template creation can use correct tuning values for resampling. This was left as-is because reversing the order degraded recognition quality for unknown reasons.

### 8. Counterexample scoring limitation

**File:** `Core/CARL/CARL/source/Recognizer.cpp`, line 418

TODO: *"This does not currently take counterexamples into account very well since they won't apply AT the site of the score."* Counterexample DTW is run against the same trimmed sequence window, but the counterexample's best match region may not overlap with the template's, reducing the penalty's effectiveness.

---

## Missing Infrastructure

### 9. Zero test coverage

No unit tests exist for `createAutoTrimmedExample`, `Definition::addExample()` round-tripping, DTW matching, or descriptor sequence construction. The `test_harness` target is a performance/integration tool using pre-baked files, not a unit test suite.

### 10. Missing emscripten binding

`createAutoTrimmedExample` is exposed in the native C API (`carl.h`, `c_api.cpp`) and Unity wrapper (`CarlExample.cs`, `CarlNative.cs`) but absent from `emscripten_bindings.cpp`. Blocker if test tooling runs in a web context (ActionStudio).

### 11. Custom descriptor type throws at runtime

**File:** `Core/CARL/CARL/source/Recognizer.cpp`, lines 183-185 (also 257-259, 305-308)

`createAutoTrimmedExample`, `analyzeRecording`, and related code paths throw `std::runtime_error{ "TODO: Not implemented" }` for `descriptor::Custom`. No compile-time guard prevents calling these on custom-type recognizers.

---

## Priority

For the adaptive definition scenario (capture user attempt → auto-trim → add as example → re-recognize):

1. **Must fix:** Timestamp boundary fallthrough (#1) — will corrupt Examples from real recordings
2. **Must understand:** Recognition pulse behavior (#3) — affects how we detect "user is performing the gesture now"
3. **Should investigate:** Delta descriptor resampling (#6) — could cause subtle trim boundary errors
4. **Should add:** Basic test coverage (#9) — at minimum, a round-trip test for auto-trim → add example → recognize
5. **Nice to have:** Emscripten binding (#10), ETS degenerate handling (#4)
