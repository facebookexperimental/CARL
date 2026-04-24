# Building Definitions from Recordings

CARL's standard workflow for defining actions requires someone to
explicitly provide `Example`s -- `Recording`s with known action
boundaries -- from which a `Definition` can be built. This works well
when the action is already known and the person providing examples can
reliably identify where the action starts and stops. But what about
scenarios where no pre-existing definition exists, and the system
itself needs to figure out what the action looks like?

Consider a guided tutorial scenario: a system instructs a user to
perform an action, then captures several `Recording`s of the user
attempting it. The system never had a `Definition` for the action --
it simply recorded whatever the user did. If the user performed
roughly the same motion each time, the system should be able to
discover the common pattern across those recordings and extract
trimmed `Example`s from them, even without knowing in advance what
to look for.

This is what `DefinitionBuilder` does. It takes multiple raw
`Recording`s of the same unknown action and produces auto-trimmed
`Example`s by finding the **common interesting sub-action** across
all recordings. These `Example`s can then be used to create a
`Definition` that encapsulates the user's own understanding of
the action.

## How It Works

`DefinitionBuilder` uses a three-step **iterative bootstrap
algorithm** to discover action boundaries without a pre-existing
`Definition` or `Recognizer`.

### Step 1: Motion-Energy Seed

For each input `Recording`, the algorithm builds a descriptor
sequence and measures **descriptor density** over time. Because
CARL's descriptor resampling is driven by motion in descriptor
space (new descriptors are emitted when the descriptor-space
distance exceeds a threshold), regions of high motion produce
denser sequences. The algorithm slides a time window across
each recording and counts descriptors per window, identifying
the region of peak activity.

The `Recording` with the most pronounced activity peak is
selected, and the densest window becomes the **seed** -- a
rough first approximation of where the action occurred.

If no significant motion is detected in any recording (as might
happen with a static gesture or pose), the algorithm selects a
half-second chunk from the middle of the recording. For static
gestures, any subsection is as representative as any other.

### Step 2: First Bootstrap Pass

The seed is used as a crude single `Example` to build a
temporary template descriptor sequence using default tuning.
Each remaining `Recording` is then auto-trimmed against this
template using Dynamic Time Warping (DTW), the same
subsequence-matching technique that `Recognizer` uses for its
own auto-trimming. This produces an initial set of trimmed
`Example`s -- one per `Recording`.

### Step 3: Refinement Pass

With multiple trimmed `Example`s now available, the algorithm
computes proper per-dimension **tuning** values from the
cross-example variance (the same tuning calculation used when
constructing a `Recognizer`). It rebuilds the template
descriptor sequences using the refined tuning, then re-auto-trims
all `Recording`s -- including the original seed -- against the
improved templates. The resulting `Example`s are returned.

## Usage

```cpp
#include <carl/DefinitionBuilder.h>
#include <carl/Definition.h>
#include <carl/Recording.h>

// Capture three recordings of the user performing an action.
// Each recording may contain leading and trailing excess.
std::vector<carl::action::Recording> recordings = captureUserRecordings();

// Discover and extract the common action from the recordings.
auto examples = carl::action::DefinitionBuilder::createExamplesFromRecordings(
    carl::action::ActionType::RightHandGesture,
    recordings,
    1.0);  // optional: expected action duration in seconds

// Build a Definition from the discovered Examples.
carl::action::Definition definition{ carl::action::ActionType::RightHandGesture };
for (auto& example : examples)
{
    definition.addExample(std::move(example));
}

// The Definition is now ready for use with a Recognizer.
```

## Parameters and Tuning

### `actionType`

The `ActionType` determines which descriptor type is used to
analyze the recordings. Choose the type that matches the kind
of action being captured:

- `LeftHandGesture` / `RightHandGesture` -- full hand gesture
  with movement (hand shape, wrist trajectory, rotation,
  translation)
- `LeftHandPose` / `RightHandPose` -- hand shape with wrist
  orientation, but without tracking movement over time
- `LeftHandShape` / `RightHandShape` -- hand shape only,
  independent of wrist orientation or movement
- `TwoHandGesture` -- gestures involving both hands
- `LeftWristTrajectory` / `RightWristTrajectory` -- wrist
  movement patterns without hand shape
- `LeftControllerGesture` / `RightControllerGesture` /
  `TwoControllerGesture` -- controller-based input

### `expectedActionDuration`

An optional hint (in seconds) for the expected length of the
action. When provided, it constrains the size of the
motion-energy window used during seeding. When set to 0 (the
default), the algorithm derives a window size from the recording
lengths (half the shortest recording).

Providing this hint improves results when you have a rough idea
of how long the action should take, especially if the recordings
contain long idle periods before or after the action.

### Number of Recordings

- **Minimum**: 2 recordings. With only 1 recording, the
  algorithm returns a single example based on the motion-energy
  seed alone (no cross-validation is possible).
- **Recommended**: 3 or more recordings. More recordings improve
  the robustness of the common-pattern discovery and the quality
  of the tuning calculation.

### Recording Length and Excess

Recordings should contain some leading and trailing time before
and after the action, but not excessively so. A good guideline
is 1-2 seconds of excess on each side. Very long recordings
(minutes of idle with a brief action) will still work but may
reduce seed quality, since the motion-energy peak becomes less
distinctive against a longer baseline.

## Limitations and Best Practices

- **Consistency matters.** The algorithm discovers the *common*
  pattern across recordings. If the user performs a substantially
  different motion each time, the extracted examples may not
  represent a coherent action. Three consistent performances are
  much more valuable than ten inconsistent ones.

- **Static gestures work but differently.** For gestures
  involving minimal motion (static hand poses), the motion-energy
  heuristic cannot distinguish the action from idle time. The
  algorithm handles this gracefully by selecting a representative
  subsection, but the quality depends on the user holding the
  pose for a consistent duration across recordings.

- **Quality validation.** After building a `Definition` from the
  extracted `Example`s, you can validate quality by creating a
  `Recognizer` and scoring the original recordings against it.
  If the recognizer scores all recordings highly, the extraction
  was successful. Low scores suggest the recordings were too
  inconsistent.

- **Coordinate conventions.** All data in CARL follows OpenXR
  conventions. Ensure that `InputSample` poses are expressed
  accordingly.
