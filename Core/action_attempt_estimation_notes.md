# Action Attempt Estimation: Design & Implementation Notes

Written April 2026. This document captures the full context of the action attempt estimation feature so work can be resumed after a hiatus.

---

## The Problem

CARL recognizes gestures by matching incoming descriptor sequences against pre-recorded templates using Dynamic Time Warping (DTW). When the match is close enough (score > 0.7), recognition fires. When it isn't, nothing happens — the system is blind to the user's intent.

This matters because the most common failure mode in gesture interaction is not "the system broke" but "the user did the gesture wrong." A user who's forgotten the exact motion, who was never properly taught, or who's physically unable to reproduce the template precisely will repeatedly attempt something *close* to the gesture without ever crossing the recognition threshold. Currently, CARL treats these attempts identically to ambient hand motion. That's a missed opportunity.

## The Motivation

The action attempt estimation feature reframes recognition: instead of asking "did the user do the gesture correctly?", CARL asks "is the user *trying* to do the gesture?" If CARL can detect repeated imperfect attempts, it can surface them to the app layer, which can then decide to adapt, guide, or accommodate.

Target scenarios include:
- A user who was shown a gesture briefly and is attempting to reproduce it from memory
- A user who has forgotten the exact motion but knows it was "something like this"
- Adaptive personalization: learning how a specific user performs a gesture over time

## The Core Insight

When a user repeatedly attempts a gesture they don't quite know:

1. The DTW **distance** (not score — distance is pre-scoring-function and retains more information) will show periodic dips below the ambient baseline, corresponding to each attempt.
2. These dips won't reach the recognition threshold, but they'll be significantly below what random ambient motion produces.
3. The dips will cluster temporally — a few seconds apart as the user resets and tries again.
4. Critically, the dips will be **similar to each other** — the user is doing the same wrong thing repeatedly.

Point 4 is the key insight that makes this tractable. The question "are these attempts similar to each other?" is much easier to answer than "are these attempts similar enough to the template to indicate intent?" because the former doesn't depend on having correctly calibrated cross-dimension normalization. You're comparing like-with-like.

## Architecture

### What CARL Owns

- **Detection**: Tracking distance history, maintaining an ambient baseline via EMA, identifying when distance drops significantly below ambient at a local minimum
- **Snapshot capture**: Caching the descriptor slice and corresponding InputSamples from each detected attempt
- **Cross-recognition**: Running DTW between cached descriptor snapshots to determine if they cluster
- **Surfacing**: Firing a signal with the clustered attempt data

### What CARL Does NOT Own

- Whether to adapt a definition from the attempts
- Whether to create a new definition
- Whether to show the user instructional feedback
- Whether to lower sensitivity
- Any policy decision about what to do with detected attempts

This is the app's responsibility. CARL provides the signal and the raw material.

## Key Design Decisions

### Distance, Not Score

The scoring function (`1 - (d / (3.16228 * sensitivity))²`) is lossy — it compresses the interesting part of the distance space. A sequence of "close but not quite" attempts might all map to scores in the 0.2-0.5 range, which is indistinguishable mush. In raw distance space, the same attempts show clear periodic dips toward a consistent value. All tracking and detection operates on distance.

### Significance Threshold Is Definition-Specific

What counts as "close" depends entirely on the gesture. A pinch gesture starting from a relaxed hand is halfway to "matched" before the user even starts. A fist-with-pinky-out gesture is miles from ambient. A single global threshold would be useless.

The threshold is either:
- **Pre-computed** at definition creation time by running unrelated motion through the recognizer (better prior, but optional and potentially stale)
- **Estimated on-the-fly** via exponential moving average of ambient distance (adapts to context, but needs warm-up)

If both are available, the pre-computed value seeds the EMA, avoiding cold-start.

### Cached Descriptors, Not Recomputed

Snapshots cache descriptor slices directly from the Provider's sequence at near-miss time — the descriptors already exist and don't need recomputation. InputSamples are also cached (from the retained sample ring buffer) but are only needed if the app wants to export/serialize the result as a Definition. The cross-recognition step operates entirely on cached descriptors.

### Use the Proximate Recognizer's Tuning

Tuning (per-dimension distance scaling) is currently effectively all 1.0 across all descriptors. The identicality thresholds baked into each descriptor's distance function do the real normalization work. When tuning values *are* set, they typically reflect broad "this action doesn't care about these dimensions" semantics (e.g., a pinch gesture that doesn't care about the other three fingers).

For cross-recognition between snapshots, using the recognizer's existing `m_tuning` is the obvious choice. It's consistent with how the recognizer evaluates the template match, and in the rare case where dimensions are muted, that muting should carry through.

### New Definitions from Snapshots, Not Expanding Existing

Adding discovered attempt examples to the existing definition would widen the tuning envelope toward the user's imperfect version, potentially causing definition drift. A separate definition built from only the discovered examples preserves the original. This is an app-level decision, but the architecture supports it by surfacing raw InputSample vectors rather than pre-built definitions.

## Implementation Summary

### Infrastructure (Phases 1 & 2 — Implemented)

**Timestamp propagation**: `extendSequence()` in `SequenceOperations.h` now requires a `std::vector<double>& timestamps` parameter. Timestamps are captured using the exact same `t` values from bisection/interpolation — not linearly interpolated after the fact. The Provider in `SessionImpl.h` maintains a parallel `m_timestamps` vector kept in sync with the descriptor sequence through trimming. `Session::Impl` exposes `getProviderTimestamps<DescriptorT>()`.

**Sample retention**: A templated `RingBuffer<T>` in `Foundation/include/carl/RingBuffer.h` provides O(1) push with ordered iteration. The Provider optionally retains InputSamples in a `RingBuffer<InputSample>`. `Session::Impl` exposes `enableSampleRetention<DescriptorT>()` and `getSamplesInRange<DescriptorT>()` (returns `std::vector<InputSample>`). The public `Session` API exposes `enableSampleRetention(capacity)` which auto-enables on Providers when they're created.

### Scoring Enrichment (Phase 3 — Implemented)

`calculateMatchDistanceIncremental` and `calculateMatchDistance` in `Recognizer.cpp` now return the full `DynamicTimeWarping::MatchResult<NumberT>` instead of a scalar. The template loop in `calculateScore` tracks the argmin template via side-effect members: `m_lastBestMatchResult`, `m_lastBestTemplateIdx`, `m_lastMinDistance`.

### ActionAttemptEstimator (Phase 4 — Implemented)

**Public types** (`ActionAttemptEstimation.h`):
- `ActionAttemptEstimationSettings`: EMA decay alpha (0.02), significance multiple (0.6), cluster window (5s), cluster min count (3), snapshot expiry (10s), optional precomputed significance threshold
- `AttemptCluster`: vector of InputSample vectors per snapshot, first/last timestamps, count

**Internal class** (`ActionAttemptEstimator.h`): Templated on `DescriptorT`.

- **EMA ambient baseline**: Updated every frame. 60-frame warmup period before detection activates (skipped if precomputed threshold is provided).
- **Local minimum detection**: 3-frame ring buffer of recent distances. A minimum at frame t-1 is confirmed when `distance[t-1] < distance[t-2] AND distance[t-1] < distance[t]`. Only minima that also clear the ambient significance threshold qualify.
- **Snapshot capture**: On significant minimum, copies the descriptor slice from the Provider's sequence (using `ImageStartIdx` and `ImageSize` from the `MatchResult`) and retrieves InputSamples from the Provider's retained sample buffer (using timestamps to find the time range).
- **Deduplication**: A "pending" snapshot is held before being committed. If a new minimum overlaps and improves on the pending snapshot, it replaces it. If an actual match occurs (score > 0.7), overlapping pending snapshots are discarded.
- **Expiry**: Committed snapshots older than `snapshotExpirySeconds` are purged.
- **Cross-DTW clustering**: When enough snapshots accumulate within the cluster window, all pairs are compared via `DynamicTimeWarping::Match` using the recognizer's distance functions and `m_tuning`. If mean pairwise distance is below `ambient × significanceMultiple`, the snapshots cluster and an `AttemptCluster` is returned.

### Wiring (Phase 5 — Implemented)

`handleSequence` in `RecognizerImpl` feeds the estimator when enabled:
1. Retrieves Provider timestamps
2. Adjusts `MatchResult.ImageStartIdx` by the trim offset (the descriptor sequence is trimmed for scoring, but timestamps/snapshots operate on the full sequence)
3. Calls `estimator->update()` with the full sequence, timestamps, session impl reference, and score
4. Calls `estimator->expireSnapshots()`
5. Calls `estimator->checkForCluster()` — if cluster detected, raises the signal

`enableActionAttemptEstimation` creates the estimator and enables sample retention on the Provider.

### Public API (Phase 6 — Implemented)

- `Recognizer::enableActionAttemptEstimation(ActionAttemptEstimationSettings)` — opt-in, creates the estimator internally
- `Recognizer::whenAttemptClusterDetectedSignal` — `Signal<AttemptCluster>`, fires when cross-recognized snapshots cluster

## Known Limitations & Future Work

### From the Auto-Trimming Health Notes (`Core/auto_trimming_health_notes.md`)

Several pre-existing issues in the auto-trimming code should be addressed:
1. **Timestamp boundary fallthrough** in `createAutoTrimmedExample` — can produce Examples with absurd timestamp bounds
2. **Recognition state machine is pulse-only** — resets after one frame regardless of score
3. **Delta descriptor interaction with resampling** — open question about correctness
4. **EgocentricTemporalSpace identity return** on degenerate geometry
5. **Zero test coverage** for all of this

### Tuning

Tuning is computed but effectively stays at ≥1.0 for all descriptors. The original purpose (per-dimension variance normalization) has been largely superseded by identicality-thresholded distance functions, multiple examples, and counterexamples. Tuning may evolve to serve a simpler purpose like dimension muting (e.g., "this gesture doesn't care about finger positions").

### Per-Dimension Analysis

The current cross-DTW clustering operates on aggregate distance. A richer signal would be the per-dimension distance breakdown at each snapshot — this would enable instructional feedback ("your hand shape is right, but your wrist trajectory is wrong"). The distance functions compute per-dimension values inside `CombinedDescriptor::AbsoluteDistance` and `DeltaDistance` but sum them before returning. Surfacing the breakdown requires either returning the vector or recomputing at snapshot time.

### Sensitivity to Non-Linear Dimension Penalties

If one descriptor dimension has a much steeper normalization curve than others, it can dominate the aggregate distance even when the user is only wrong in that dimension. This could cause some users' attempts to register as significant near-misses while others' don't, even if both users forgot the same number of gesture elements. The cross-recognition approach (comparing attempts to each other rather than to templates) partially mitigates this, but the effect should be monitored.

### C API & Bindings

The action attempt estimation feature is not yet exposed through the C API, Unity wrapper, or emscripten bindings. This is needed before it can be used in production contexts.

### Cluster Threshold Calibration

The current cluster threshold (`ambient × significanceMultiple`) is a reasonable starting point but has not been empirically validated. Real-world testing with actual user data will likely require tuning the defaults.

## File Map

| File | Role |
|------|------|
| `Core/CARL/Utilities/include/carl/descriptor/SequenceOperations.h` | Timestamp propagation in `extendSequence` |
| `Core/CARL/Foundation/include/carl/RingBuffer.h` | General-purpose ring buffer utility |
| `Core/CARL/CARL/source/SessionImpl.h` | Provider timestamp storage, sample retention, Session::Impl accessors |
| `Core/CARL/CARL/include/carl/Session.h` | Public `enableSampleRetention` API |
| `Core/CARL/CARL/source/Session.cpp` | `enableSampleRetention` forwarding |
| `Core/CARL/CARL/include/carl/ActionAttemptEstimation.h` | Public types: Settings, AttemptCluster |
| `Core/CARL/CARL/source/ActionAttemptEstimator.h` | Internal estimator: EMA, detection, snapshots, cross-DTW |
| `Core/CARL/CARL/include/carl/Recognizer.h` | Public API: `enableActionAttemptEstimation`, cluster signal |
| `Core/CARL/CARL/source/Recognizer.cpp` | Scoring enrichment, estimator wiring, API forwarding |
| `Core/auto_trimming_health_notes.md` | Pre-existing issues in auto-trimming code |
