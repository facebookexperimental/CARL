# Tutorial: Recognizing Actions

This tutorial shows how to use a `Definition` for live action recognition — creating sessions, feeding input, and responding to recognition events.

**Prerequisites:** Familiarity with [Concepts](../Concepts.md) and [Defining Actions](DefiningActions.md).

## Step 1: Creating a Session

A `Session` is the runtime context that receives live input and hosts recognizers:

```cpp
#include <carl/Carl.h>

using namespace carl;
using namespace carl::action;

// Multi-threaded (default) — recognition runs on a background thread
Session session{};

// Single-threaded — everything runs on the calling thread
Session session{/* singleThreaded */ true};
```

**Multi-threaded mode** (default) is recommended for production. Recognition processing happens on a background thread, and results are dispatched to the main thread via `tickCallbacks()`.

**Single-threaded mode** is useful for debugging or environments where threading is not available. All processing happens inline during `addInput()`.

### Logging

Enable diagnostic logging to help debug recognition issues:

```cpp
session.setLogger([](std::string msg) {
    std::cout << "[CARL] " << msg << std::endl;
});
```

## Step 2: Feeding Input

Feed `InputSample`s to the session every frame:

```cpp
void onFrame(double timestamp, /* ... input data ... */) {
    InputSample sample{};
    sample.Timestamp = timestamp;
    sample.HmdPose = getHmdPose();
    sample.RightWristPose = getRightWristPose();
    sample.RightHandJointPoses = getRightHandJoints();

    session.addInput(sample);
}
```

Input should be provided at your application's frame rate. CARL does not require a specific frequency, but faster input generally produces smoother scores.

## Step 3: Creating a Recognizer

A `Recognizer` scores live input against a `Definition`. Create one for each action you want to recognize:

```cpp
// Load a saved definition
auto loaded = carl::utilities::TryDeserializeFromFile<Definition>("wave.carl");
Definition& definition = loaded.value();

// Create the recognizer
Recognizer recognizer{session, definition};
```

You can create multiple recognizers on the same session — each independently evaluates its own definition against the shared input stream.

## Step 4: Reading Results

There are two ways to get recognition results: polling and event-based.

### Option A: Polling

Check the current score each frame:

```cpp
void onFrame() {
    double score = recognizer.currentScore();

    if (score >= 1.0) {
        // Action is being recognized
    }
}
```

`currentScore()` is thread-safe and can be called from any thread at any time.

### Option B: Event-Based

Subscribe to `whenRecognitionChangedSignal` for discrete start/stop notifications:

```cpp
auto ticket = recognizer.whenRecognitionChangedSignal.addHandler(
    [](bool recognized) {
        if (recognized) {
            std::cout << "Action started!" << std::endl;
        } else {
            std::cout << "Action ended." << std::endl;
        }
    });
```

The signal fires `true` when the score crosses above the recognition threshold and `false` when it drops below. The handler is called from `tickCallbacks()`, so it runs on the thread that calls `tickCallbacks()`.

**Important:** The returned `ticket` controls the subscription lifetime. Store it for as long as you want to receive events. When the ticket is destroyed, the handler is automatically unsubscribed.

## Step 5: Ticking Callbacks

In multi-threaded mode, recognition happens on a background thread. Call `tickCallbacks()` on your main thread to dispatch signal handlers:

```cpp
arcana::cancellation token{};

void onFrame() {
    session.addInput(sample);
    session.tickCallbacks(token);
    // Signal handlers fire here, on this thread
}
```

In single-threaded mode, callbacks are dispatched inline during `addInput()`, so `tickCallbacks()` is optional but harmless.

## Step 6: Canonical Recording

Get the recognizer's internal view of the input for debugging or visualization:

```cpp
RecordingInspector inspector = recognizer.getCanonicalRecordingInspector();
double start = inspector.startTimestamp();
double end = inspector.endTimestamp();

// Inspect the input at a specific time
const InputSample& sample = inspector.inspect(someTimestamp);
```

This can be useful for understanding what the recognizer "sees" and for building visualization tools.

## Complete Working Example

```cpp
#include <carl/Carl.h>
#include <carl/utilities/FileSerialization.h>
#include <iostream>

using namespace carl;
using namespace carl::action;

int main() {
    // 1. Create session
    Session session{/* singleThreaded */ true};
    session.setLogger([](std::string msg) {
        std::cout << "[CARL] " << msg << std::endl;
    });

    // 2. Load a definition
    auto loaded = carl::utilities::TryDeserializeFromFile<Definition>("wave.carl");
    if (!loaded.has_value()) {
        std::cerr << "Failed to load definition" << std::endl;
        return 1;
    }
    Definition& definition = loaded.value();

    // 3. Create a recognizer
    Recognizer recognizer{session, definition};

    // 4. Subscribe to events
    auto ticket = recognizer.whenRecognitionChangedSignal.addHandler(
        [](bool recognized) {
            std::cout << (recognized ? "Wave detected!" : "Wave ended.") << std::endl;
        });

    // 5. Main loop — feed input and tick
    arcana::cancellation token{};
    while (/* running */) {
        InputSample sample{};
        sample.Timestamp = getCurrentTimeInSeconds();
        // Populate sample fields from XR hardware...

        session.addInput(sample);
        session.tickCallbacks(token);

        // Optionally poll the score
        double score = recognizer.currentScore();
        std::cout << "Score: " << score << std::endl;
    }

    return 0;
}
```

## Multiple Recognizers

You can recognize multiple actions simultaneously:

```cpp
auto waveDef = carl::utilities::TryDeserializeFromFile<Definition>("wave.carl").value();
auto snapDef = carl::utilities::TryDeserializeFromFile<Definition>("snap.carl").value();

Recognizer waveRecognizer{session, waveDef};
Recognizer snapRecognizer{session, snapDef};

auto waveTicket = waveRecognizer.whenRecognitionChangedSignal.addHandler(
    [](bool r) { if (r) std::cout << "Wave!" << std::endl; });
auto snapTicket = snapRecognizer.whenRecognitionChangedSignal.addHandler(
    [](bool r) { if (r) std::cout << "Snap!" << std::endl; });
```

All recognizers share the same input stream from the session.

## Next Steps

- [Tuning and Debugging](TuningAndDebugging.md) — improve recognition quality
- [C++ API Reference](../CppApiReference.md) — full API details
- [Serialization](../Serialization.md) — save and load definitions
