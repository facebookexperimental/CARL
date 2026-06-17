# C++ API Reference

Complete reference for the CARL C++ API. All types are in the `carl` or `carl::action` namespace unless otherwise noted.

Include `<carl/Carl.h>` for all types, or include individual headers as listed below.

## ActionType

**Header:** `<carl/ActionType.h>`
**Namespace:** `carl::action`

```cpp
enum class ActionType : uint64_t {
    LeftHandPose            = 0,
    LeftHandGesture         = 1,
    RightHandPose           = 2,
    RightHandGesture        = 3,
    TwoHandGesture          = 4,
    LeftControllerGesture   = 5,
    RightControllerGesture  = 6,
    TwoControllerGesture    = 7,
    LeftWristTrajectory     = 8,
    RightWristTrajectory    = 9,
    LeftHandShape           = 10,
    RightHandShape          = 11,
    Custom                  = ~uint64_t{0x0}
};
```

Determines which `InputSample` fields are used and how they are compared during recognition. See [Action Types](ActionTypes.md) for detailed descriptions of each value.

## InputSample

**Header:** `<carl/InputSample.h>`
**Namespace:** `carl`

A snapshot of all XR input state at a single instant. All 3D poses must be in the right-handed Y-up (OpenXR) coordinate system.

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `Timestamp` | `double` | Time in seconds (required) |
| `HmdPose` | `std::optional<TransformT>` | Head-mounted display 6DOF pose |
| `LeftWristPose` | `std::optional<TransformT>` | Left wrist 6DOF pose |
| `RightWristPose` | `std::optional<TransformT>` | Right wrist 6DOF pose |
| `LeftHandJointPoses` | `std::optional<std::array<TransformT, 26>>` | 26 left-hand joint poses |
| `RightHandJointPoses` | `std::optional<std::array<TransformT, 26>>` | 26 right-hand joint poses |
| `LeftControllerInput` | `std::optional<std::array<NumberT, 7>>` | 7 left controller input channels |
| `RightControllerInput` | `std::optional<std::array<NumberT, 7>>` | 7 right controller input channels |

### Nested Enum: Joint

```cpp
enum class Joint : uint64_t {
    XR_HAND_JOINT_PALM_EXT                  = 0,
    XR_HAND_JOINT_WRIST_EXT                 = 1,
    XR_HAND_JOINT_THUMB_METACARPAL_EXT      = 2,
    XR_HAND_JOINT_THUMB_PROXIMAL_EXT        = 3,
    XR_HAND_JOINT_THUMB_DISTAL_EXT          = 4,
    XR_HAND_JOINT_THUMB_TIP_EXT             = 5,
    XR_HAND_JOINT_INDEX_METACARPAL_EXT      = 6,
    XR_HAND_JOINT_INDEX_PROXIMAL_EXT        = 7,
    XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT    = 8,
    XR_HAND_JOINT_INDEX_DISTAL_EXT          = 9,
    XR_HAND_JOINT_INDEX_TIP_EXT             = 10,
    XR_HAND_JOINT_MIDDLE_METACARPAL_EXT     = 11,
    XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT       = 12,
    XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT   = 13,
    XR_HAND_JOINT_MIDDLE_DISTAL_EXT         = 14,
    XR_HAND_JOINT_MIDDLE_TIP_EXT            = 15,
    XR_HAND_JOINT_RING_METACARPAL_EXT       = 16,
    XR_HAND_JOINT_RING_PROXIMAL_EXT         = 17,
    XR_HAND_JOINT_RING_INTERMEDIATE_EXT     = 18,
    XR_HAND_JOINT_RING_DISTAL_EXT           = 19,
    XR_HAND_JOINT_RING_TIP_EXT             = 20,
    XR_HAND_JOINT_LITTLE_METACARPAL_EXT     = 21,
    XR_HAND_JOINT_LITTLE_PROXIMAL_EXT       = 22,
    XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT   = 23,
    XR_HAND_JOINT_LITTLE_DISTAL_EXT         = 24,
    XR_HAND_JOINT_LITTLE_TIP_EXT            = 25,
    COUNT                                    = 26
};
```

Follows the [OpenXR XR_EXT_hand_tracking](https://registry.khronos.org/OpenXR/specs/1.0/html/xrspec.html#XR_EXT_hand_tracking) joint layout.

### Nested Enum: ControllerInput

```cpp
enum class ControllerInput : uint64_t {
    XR_KHR_generic_controller_PRIMARY_CLICK      = 0,
    XR_KHR_generic_controller_SECONDARY_CLICK    = 1,
    XR_KHR_generic_controller_THUMBSTICK_X       = 2,
    XR_KHR_generic_controller_THUMBSTICK_Y       = 3,
    XR_KHR_generic_controller_THUMBSTICK_CLICK   = 4,
    XR_KHR_generic_controller_SQUEEZE_VALUE      = 5,
    XR_KHR_generic_controller_TRIGGER_VALUE      = 6,
    COUNT                                         = 7
};
```

### Methods

| Method | Description |
|--------|-------------|
| `InputSample()` | Default constructor; all optional fields are empty |
| `InputSample(Deserialization&)` | Deserializing constructor |
| `void serialize(Serialization&) const` | Serialize to binary |
| `static InputSample Lerp(const InputSample& a, const InputSample& b, double t)` | Linearly interpolate between two samples at parameter `t` ∈ [0, 1] |

## InProgressRecording

**Header:** `<carl/Recording.h>`
**Namespace:** `carl::action`

A mutable, rolling buffer of `InputSample`s. Move it into a `Recording` to finalize.

### Constructors

| Constructor | Description |
|-------------|-------------|
| `InProgressRecording()` | Unlimited buffer — keeps all samples |
| `InProgressRecording(size_t maxSeconds)` | Rolling buffer — discards samples older than `maxSeconds` |

Copy is deleted; move is allowed.

### Methods

| Method | Description |
|--------|-------------|
| `void addSample(InputSample sample)` | Appends a sample. In rolling mode, drops samples older than `maxSeconds` from the current timestamp. |

## Recording

**Header:** `<carl/Recording.h>`
**Namespace:** `carl::action`

An immutable, finalized sequence of `InputSample`s.

### Constructors

| Constructor | Description |
|-------------|-------------|
| `Recording(InProgressRecording)` | Consumes and finalizes an `InProgressRecording` |
| `Recording(Deserialization&)` | Deserializing constructor |

Copy is allowed; move is allowed. Assignment is deleted.

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `serialize(Serialization&) const` | `void` | Serialize to binary |
| `getInspector() const` | `RecordingInspector` | Returns a non-owning inspector view. The `Recording` must outlive the inspector. |
| `getSamples() const` | `gsl::span<const InputSample>` | Direct access to the raw sample array |

## RecordingInspector

**Header:** `<carl/Recording.h>`
**Namespace:** `carl::action`

A non-owning view for efficient timestamp-based scrubbing through a recording. The underlying data source must outlive the inspector.

### Constructor

| Constructor | Description |
|-------------|-------------|
| `RecordingInspector(gsl::span<const InputSample>)` | Creates an inspector over the given samples |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `inspect(double timestamp)` | `const InputSample&` | Returns the sample closest to `timestamp` |
| `startTimestamp() const` | `double` | Timestamp of the first sample |
| `endTimestamp() const` | `double` | Timestamp of the last sample |

## Example

**Header:** `<carl/Example.h>`
**Namespace:** `carl::action`

A `Recording` annotated with start and end timestamps marking when an action occurred.

### Constructors

| Constructor | Description |
|-------------|-------------|
| `Example(Recording recording, double startTimestamp, double endTimestamp)` | Creates an example from a recording with action boundaries |
| `Example(Deserialization&)` | Deserializing constructor |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `serialize(Serialization&) const` | `void` | Serialize to binary |
| `getRecording() const` | `const Recording&` | Access the underlying recording |
| `getStartTimestamp() const` | `double` | Get the action start time |
| `setStartTimestamp(double)` | `void` | Update the action start time |
| `getEndTimestamp() const` | `double` | Get the action end time |
| `setEndTimestamp(double)` | `void` | Update the action end time |

## Definition

**Header:** `<carl/Definition.h>`
**Namespace:** `carl::action`

A complete description of a recognizable action: one or more examples, optional counterexamples, an `ActionType`, and a sensitivity value.

### Constructors

| Constructor | Description |
|-------------|-------------|
| `Definition(ActionType actionType)` | Creates an empty definition for the given action type |
| `Definition(Deserialization&)` | Deserializing constructor |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `serialize(Serialization&) const` | `void` | Serialize to binary |
| `addExample(Example)` | `void` | Add a positive example |
| `addCounterexample(Example)` | `void` | Add a negative example for disambiguation |
| `getDescriptorType() const` | `ActionType` | Get the action type |
| `getExamples() const` | `gsl::span<const Example>` | Access all positive examples |
| `getCounterexamples() const` | `gsl::span<const Example>` | Access all counterexamples |

### Public Data Members

| Member | Type | Default | Description |
|--------|------|---------|-------------|
| `DefaultSensitivity` | `double` | `1.0` | Default recognition sensitivity. Higher values make the action easier to trigger. |

## Session

**Header:** `<carl/Session.h>`
**Namespace:** `carl`

The runtime context for action recognition. Receives live input and hosts recognizers.

### Type Aliases

```cpp
using SchedulerT = stdext::inplace_function<void(stdext::inplace_function<void(), 128>&&), 128>;
```

### Constructor

| Constructor | Description |
|-------------|-------------|
| `Session(bool singleThreaded = false)` | Creates a session. Pass `true` for single-threaded mode. |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `addInput(const InputSample&)` | `void` | Feed a live input sample into the session |
| `tickCallbacks(arcana::cancellation& token)` | `void` | Dispatch queued signal handlers on the calling thread |
| `setLogger(std::function<void(std::string)>)` | `void` | Set a logging callback for diagnostics |
| `log(std::string message)` | `void` | Emit a log message through the logger |
| `callbackScheduler()` | `SchedulerT&` | Access the callback scheduler (advanced) |
| `processingScheduler()` | `SchedulerT&` | Access the processing scheduler (advanced) |
| `enableCustomActionType<T>(TemplatedCustomActionTypeOperations<T>)` | `ContractId<>::IdT` | Register a custom action type; returns an ID for use with `Recognizer` |

### TemplatedCustomActionTypeOperations\<T\>

Template struct for registering custom action types:

```cpp
template<typename T>
struct TemplatedCustomActionTypeOperations {
    // Create a descriptor from current and prior InputSamples.
    // Return std::nullopt if required data is unavailable.
    std::function<std::optional<T>(const InputSample& current,
                                    const InputSample& prior)> TryCreate;

    // Distance between two descriptors (a, b) considering their origins (a0, b0)
    // and per-dimension tuning weights.
    std::function<NumberT(const T& a, const T& a0,
                           const T& b, const T& b0,
                           gsl::span<const NumberT> tuning)> Distance;

    // Linear interpolation between two descriptors.
    std::function<T(const T& a, const T& b, NumberT t)> Lerp;

    // Calculate per-dimension tuning weights from training examples.
    std::function<std::array<NumberT, 32>(
        gsl::span<const action::Example> examples)> CalculateTuning;
};
```

## Recognizer

**Header:** `<carl/Recognizer.h>`
**Namespace:** `carl::action`

Observes a `Session` and produces recognition scores for a `Definition`.

### Constructors

| Constructor | Description |
|-------------|-------------|
| `Recognizer(Session&, const Definition&)` | Creates a recognizer for a built-in action type |
| `Recognizer(Session&, const Definition&, ContractId<>::IdT)` | Creates a recognizer for a custom action type, using the ID returned by `enableCustomActionType()` |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `currentScore() const` | `double` | Current recognition confidence. ≥ 1.0 = recognized, ≈ 0.0 = not recognized. Thread-safe. |
| `getSensitivity() const` | `double` | Get the current sensitivity |
| `setSensitivity(double)` | `void` | Override the definition's default sensitivity |
| `getCanonicalRecordingInspector() const` | `RecordingInspector` | Inspector for the recognizer's internal canonical representation |
| `createAutoTrimmedExample(const Recording&) const` | `Example` | Automatically determine optimal start/end timestamps for a new example |
| `analyzeRecording(const Recording&, std::ostream&) const` | `void` | Output diagnostic analysis of how a recording matches the definition |

### Public Data Members

| Member | Type | Description |
|--------|------|-------------|
| `whenRecognitionChangedSignal` | `Signal<bool>` | Fires `true` when action starts being recognized, `false` when it stops |

## Signal\<...ArgsT\>

**Header:** `<carl/Signaling.h>`
**Namespace:** `carl`

A lightweight publish/subscribe signal.

### Type Aliases

```cpp
using HandlerT = stdext::inplace_function<void(ArgsT&...), 128>;
using TicketT  = typename arcana::weak_table<HandlerT>::ticket;
```

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `addHandler(CallableT&& callable)` | `TicketT` | Subscribe a handler. Returns a ticket that controls subscription lifetime — destroying the ticket unsubscribes the handler. |

## Serialization / Deserialization

**Header:** `<carl/Serialization.h>`
**Namespace:** `carl`

Low-level binary serialization. Most users should prefer the file utilities below.

### Serialization

```cpp
Serialization(std::vector<uint8_t>& bytes); // Clears bytes and writes into it
```

Provides `operator<<` for POD types, `std::vector`, `std::array`, `std::optional`, `std::string`, and CARL math types.

### Deserialization

```cpp
Deserialization(const uint8_t* bytes, size_t length);
Deserialization(gsl::span<const uint8_t> bytes);
```

Provides `operator>>` for the same types. Also: `size_t remainingBytes() const`.

## File Serialization Utilities

**Header:** `<carl/utilities/FileSerialization.h>`
**Namespace:** `carl::utilities`

High-level functions for serializing `Example` and `Definition` to/from bytes and files.

| Function | Description |
|----------|-------------|
| `Serialize<T>(const T&)` | Serialize to `std::vector<uint8_t>` |
| `SerializeToFile<T>(const T&, std::filesystem::path)` | Serialize to a file |
| `TryDeserialize<T>(gsl::span<const uint8_t>)` | Deserialize from bytes; returns `std::optional<T>` |
| `TryDeserializeFromFile<T>(std::filesystem::path)` | Deserialize from a file; returns `std::optional<T>` |
| `TryDeserializeLegacy<T>(gsl::span<const uint8_t>)` | Deserialize from legacy format |
| `TryDeserializeLegacyFile<T>(std::filesystem::path)` | Deserialize from legacy format file |

Explicit specializations are provided for `action::Example` and `action::Definition`.

## Type Aliases

**Header:** `<carl/Types.h>`
**Namespace:** `carl`

| Alias | Underlying Type |
|-------|-----------------|
| `NumberT` | `float` |
| `VectorT` | `Eigen::Vector<float, 3>` |
| `QuaternionT` | `Eigen::Quaternion<float>` |
| `AngleAxisT` | `Eigen::AngleAxis<float>` |
| `TransformT` | `Eigen::Transform<float, 3, Affine>` |

### Math Utilities (namespace `carl::math`)

| Function | Description |
|----------|-------------|
| `Lerp(float a, float b, NumberT t)` | Linear interpolation between scalars |
| `Lerp(const TransformT& a, const TransformT& b, NumberT t)` | Linear interpolation between transforms |
| `LookTransform(const VectorT& forward, const VectorT& up, const VectorT& position)` | Construct a look-at transform |
