# C API Reference

Complete reference for the CARL C API, designed for FFI from Unity (C#), WebAssembly, and other languages.

**Header:** `<carl.h>`

All functions use the `carl_` prefix. Handles are `uint64_t` values representing pointers to C++ objects.

## Conventions

### Handle Types

All CARL objects are represented as `uint64_t` handles. Functions that create objects return handles; functions that consume or inspect objects take handles as parameters.

### Ownership and Lifetime

- **Create/start/finish/deserialize** functions allocate objects — you own the result and must call the corresponding `dispose` function when done.
- **Dispose** functions free the object. Using a handle after disposal is undefined behavior.
- **Get** functions (e.g., `carl_getRecording`) return new handles that you must also dispose.

### Byte Buffer Pattern

Several functions return serialized data as byte buffers. These use a two-call pattern:

1. Call with `destination = NULL` and `size = 0` to query the buffer size (returned as `uint64_t`).
2. Allocate a buffer of that size.
3. Call again with the buffer pointer and size to copy the data. The native buffer is freed after this call.

```c
// Example: serialize a recording
uint64_t bytesHandle = carl_serializeRecording(recordingPtr);

// Query size
uint64_t size = carl_getBytes(bytesHandle, NULL, 0);

// Copy and free
uint8_t* buffer = malloc(size);
carl_getBytes(bytesHandle, buffer, size);
```

### InputSample Struct

The C API uses `carl_InputSample`, a struct matching the layout of the C++ `InputSample`:

```c
typedef struct {
    double Timestamp;
    OptionalTransform HmdPose;
    OptionalTransform LeftWristPose;
    OptionalTransform RightWristPose;
    OptionalTransform LeftHandJointPoses[26];
    OptionalTransform RightHandJointPoses[26];
    OptionalControllerState LeftControllerState;
    OptionalControllerState RightControllerState;
} carl_InputSample;
```

Where:

```c
typedef struct {
    bool Valid;
    struct { double X, Y, Z; } Position;
    struct { double W, X, Y, Z; } Orientation;
} OptionalTransform;

typedef struct {
    bool Valid;
    double PrimaryClick;
    double SecondaryClick;
    double ThumbstickX;
    double ThumbstickY;
    double ThumbstickClick;
    double SqueezeValue;
    double TriggerValue;
} OptionalControllerState;
```

## Byte Buffer Functions

### carl_getBytes

```c
uint64_t carl_getBytes(uint64_t bytesPtr, uint8_t* destination, uint64_t size);
```

Copies serialized bytes to `destination` and frees the native buffer. If `destination` is NULL, returns the size without copying. See [Byte Buffer Pattern](#byte-buffer-pattern).

## Recording Functions

### carl_startRecording

```c
uint64_t carl_startRecording(uint64_t maxSeconds);
```

Creates a new `InProgressRecording` with a rolling window of `maxSeconds`. Returns an `InProgressRecording` handle.

### carl_recordObjectInputSample

```c
void carl_recordObjectInputSample(uint64_t inProgressRecordingPtr, carl_InputSample* sample);
```

Adds an `InputSample` to an in-progress recording using the struct directly.

### carl_recordInputSample

```c
void carl_recordInputSample(uint64_t inProgressRecordingPtr, uint8_t* bytes, uint64_t size);
```

Adds an `InputSample` to an in-progress recording from serialized bytes.

### carl_finishRecording

```c
uint64_t carl_finishRecording(uint64_t inProgressRecordingPtr);
```

Finalizes an `InProgressRecording` into an immutable `Recording`. The `InProgressRecording` handle is consumed and should not be used after this call. Returns a `Recording` handle.

### carl_serializeRecording

```c
uint64_t carl_serializeRecording(uint64_t recordingPtr);
```

Serializes a `Recording` to bytes. Returns a byte buffer handle — use `carl_getBytes` to retrieve the data.

### carl_deserializeRecording

```c
uint64_t carl_deserializeRecording(uint8_t* bytes, uint64_t size);
```

Deserializes a `Recording` from bytes. Returns a `Recording` handle.

### carl_disposeRecording

```c
void carl_disposeRecording(uint64_t recordingPtr);
```

Frees a `Recording`.

### carl_getRecordingInspector

```c
uint64_t carl_getRecordingInspector(uint64_t recordingPtr);
```

Creates a `RecordingInspector` for a recording. Returns a `RecordingInspector` handle. The recording must outlive the inspector.

## RecordingInspector Functions

### carl_inspect

```c
uint64_t carl_inspect(uint64_t recordingInspectorPtr, double timestamp);
```

Returns a byte buffer handle containing the serialized `InputSample` closest to `timestamp`. Use `carl_getBytes` to retrieve the data.

### carl_getStartTimestamp

```c
double carl_getStartTimestamp(uint64_t recordingInspectorPtr);
```

Returns the timestamp of the first sample in the recording.

### carl_getEndTimestamp

```c
double carl_getEndTimestamp(uint64_t recordingInspectorPtr);
```

Returns the timestamp of the last sample in the recording.

### carl_disposeRecordingInspector

```c
void carl_disposeRecordingInspector(uint64_t recordingInspectorPtr);
```

Frees a `RecordingInspector`.

## Example Functions

### carl_createExample

```c
uint64_t carl_createExample(uint64_t recordingPtr, double startTimestamp, double endTimestamp);
```

Creates an `Example` from a recording with the specified action boundaries. Returns an `Example` handle.

### carl_createAutoTrimmedExample

```c
uint64_t carl_createAutoTrimmedExample(uint64_t recognizerPtr, uint64_t recordingPtr);
```

Creates an auto-trimmed `Example` using a recognizer to determine optimal start/end timestamps. Returns an `Example` handle.

### carl_getRecording

```c
uint64_t carl_getRecording(uint64_t examplePtr);
```

Gets the recording from an example. Returns a new `Recording` handle that you must dispose.

### carl_getExampleStartTimestamp

```c
double carl_getExampleStartTimestamp(uint64_t examplePtr);
```

Returns the start timestamp of the example.

### carl_getExampleEndTimestamp

```c
double carl_getExampleEndTimestamp(uint64_t examplePtr);
```

Returns the end timestamp of the example.

### carl_disposeExample

```c
void carl_disposeExample(uint64_t examplePtr);
```

Frees an `Example`.

### carl_loadExampleFromFile

```c
uint64_t carl_loadExampleFromFile(const char* path);
```

Loads an `Example` from a `.carl` file. Returns an `Example` handle, or 0 on failure.

### carl_saveExampleToFile

```c
void carl_saveExampleToFile(uint64_t examplePtr, const char* path);
```

Saves an `Example` to a `.carl` file.

## Definition Functions

### carl_createDefinition

```c
uint64_t carl_createDefinition(uint64_t descriptorType);
```

Creates an empty `Definition` for the given `ActionType` value (as `uint64_t`). Returns a `Definition` handle.

### carl_addExample

```c
void carl_addExample(uint64_t definitionPtr, uint64_t recordingPtr,
                     double startTimestamp, double endTimestamp);
```

Adds a positive example to a definition, creating an `Example` internally from the recording and timestamps.

### carl_addCounterexample

```c
void carl_addCounterexample(uint64_t definitionPtr, uint64_t recordingPtr,
                            double startTimestamp, double endTimestamp);
```

Adds a counterexample to a definition for disambiguation.

### carl_getDefaultSensitivity

```c
double carl_getDefaultSensitivity(uint64_t definitionPtr);
```

Returns the definition's default sensitivity value.

### carl_setDefaultSensitivity

```c
void carl_setDefaultSensitivity(uint64_t definitionPtr, double defaultSensitivity);
```

Sets the definition's default sensitivity value.

### carl_serializeDefinition

```c
uint64_t carl_serializeDefinition(uint64_t definitionPtr);
```

Serializes a `Definition` to bytes. Returns a byte buffer handle.

### carl_deserializeDefinition

```c
uint64_t carl_deserializeDefinition(uint8_t* bytes, uint64_t size);
```

Deserializes a `Definition` from bytes. Returns a `Definition` handle.

### carl_loadDefinitionFromFile

```c
uint64_t carl_loadDefinitionFromFile(const char* path);
```

Loads a `Definition` from a `.carl` file. Returns a `Definition` handle, or 0 on failure.

### carl_saveDefinitionToFile

```c
void carl_saveDefinitionToFile(uint64_t definitionPtr, const char* path);
```

Saves a `Definition` to a `.carl` file.

### carl_getExamplesCount

```c
uint64_t carl_getExamplesCount(uint64_t definitionPtr);
```

Returns the number of positive examples in the definition.

### carl_getExampleAtIdx

```c
uint64_t carl_getExampleAtIdx(uint64_t definitionPtr, uint64_t idx);
```

Returns an `Example` handle for the example at the given index. You must dispose the returned handle.

### carl_getCounterexamplesCount

```c
uint64_t carl_getCounterexamplesCount(uint64_t definitionPtr);
```

Returns the number of counterexamples in the definition.

### carl_getCounterexampleAtIdx

```c
uint64_t carl_getCounterexampleAtIdx(uint64_t definitionPtr, uint64_t idx);
```

Returns an `Example` handle for the counterexample at the given index. You must dispose the returned handle.

### carl_disposeDefinition

```c
void carl_disposeDefinition(uint64_t definitionPtr);
```

Frees a `Definition`.

## Session Functions

### carl_createSession

```c
uint64_t carl_createSession();
```

Creates a multi-threaded `Session`. Recognition processing runs on a background thread. Returns a `Session` handle.

### carl_createSingleThreadedSession

```c
uint64_t carl_createSingleThreadedSession();
```

Creates a single-threaded `Session`. All processing runs inline. Returns a `Session` handle.

### carl_setSessionLogger

```c
void carl_setSessionLogger(uint64_t sessionPtr, void (__cdecl* callback)(const char*));
```

Sets a logging callback for diagnostics. The callback receives null-terminated log messages. Pass `NULL` to disable logging.

### carl_tickCallbacks

```c
void carl_tickCallbacks(uint64_t sessionPtr);
```

Dispatches queued callbacks on the calling thread. Call this every frame in your main loop.

### carl_addSerializedInputSample

```c
void carl_addSerializedInputSample(uint64_t sessionPtr, uint8_t* bytes, uint64_t size);
```

Feeds a serialized `InputSample` into the session.

### carl_addInputSample

```c
void carl_addInputSample(uint64_t sessionPtr, carl_InputSample* sample);
```

Feeds an `InputSample` struct into the session.

### carl_disposeSession

```c
void carl_disposeSession(uint64_t sessionPtr);
```

Frees a `Session` and all its associated resources.

## Recognizer Functions

### carl_createRecognizerAsync

```c
void carl_createRecognizerAsync(uint64_t sessionPtr, uint64_t definitionPtr,
                                uint64_t requestId,
                                void (__cdecl* callback)(uint64_t requestId,
                                                         uint64_t recognizerPtr));
```

Asynchronously creates a `Recognizer`. The callback is called (potentially from a background thread) when creation completes, receiving the `requestId` and a `Recognizer` handle. Use `tickCallbacks` to marshal the callback to the main thread.

### carl_getCurrentScore

```c
double carl_getCurrentScore(uint64_t recognizerPtr);
```

Returns the current recognition confidence score. Thread-safe — can be called from any thread at any time.

- **≥ 1.0** — action is recognized
- **≈ 0.0** — action is not recognized

### carl_setSensitivity

```c
void carl_setSensitivity(uint64_t recognizerPtr, double sensitivity);
```

Overrides the definition's default sensitivity for this recognizer.

### carl_getCanonicalRecordingInspector

```c
uint64_t carl_getCanonicalRecordingInspector(uint64_t recognizerPtr);
```

Returns a `RecordingInspector` handle for the recognizer's internal canonical representation. You must dispose the returned handle.

### carl_disposeRecognizer

```c
void carl_disposeRecognizer(uint64_t sessionPtr, uint64_t recognizerPtr);
```

Frees a `Recognizer`. Requires the session handle because recognizer disposal is coordinated with the session's processing thread.

## Function Summary

| Group | Function | Action |
|-------|----------|--------|
| Bytes | `carl_getBytes` | Copy/query byte buffer |
| Recording | `carl_startRecording` | Create InProgressRecording |
| Recording | `carl_recordObjectInputSample` | Add sample (struct) |
| Recording | `carl_recordInputSample` | Add sample (bytes) |
| Recording | `carl_finishRecording` | Finalize recording |
| Recording | `carl_serializeRecording` | Serialize |
| Recording | `carl_deserializeRecording` | Deserialize |
| Recording | `carl_disposeRecording` | Free |
| Recording | `carl_getRecordingInspector` | Create inspector |
| Inspector | `carl_inspect` | Get sample at timestamp |
| Inspector | `carl_getStartTimestamp` | First timestamp |
| Inspector | `carl_getEndTimestamp` | Last timestamp |
| Inspector | `carl_disposeRecordingInspector` | Free |
| Example | `carl_createExample` | Create from recording |
| Example | `carl_createAutoTrimmedExample` | Auto-trim via recognizer |
| Example | `carl_getRecording` | Get recording |
| Example | `carl_getExampleStartTimestamp` | Get start time |
| Example | `carl_getExampleEndTimestamp` | Get end time |
| Example | `carl_disposeExample` | Free |
| Example | `carl_loadExampleFromFile` | Load from file |
| Example | `carl_saveExampleToFile` | Save to file |
| Definition | `carl_createDefinition` | Create empty |
| Definition | `carl_addExample` | Add positive example |
| Definition | `carl_addCounterexample` | Add negative example |
| Definition | `carl_getDefaultSensitivity` | Get sensitivity |
| Definition | `carl_setDefaultSensitivity` | Set sensitivity |
| Definition | `carl_serializeDefinition` | Serialize |
| Definition | `carl_deserializeDefinition` | Deserialize |
| Definition | `carl_loadDefinitionFromFile` | Load from file |
| Definition | `carl_saveDefinitionToFile` | Save to file |
| Definition | `carl_getExamplesCount` | Count examples |
| Definition | `carl_getExampleAtIdx` | Get example by index |
| Definition | `carl_getCounterexamplesCount` | Count counterexamples |
| Definition | `carl_getCounterexampleAtIdx` | Get counterexample by index |
| Definition | `carl_disposeDefinition` | Free |
| Session | `carl_createSession` | Create multi-threaded |
| Session | `carl_createSingleThreadedSession` | Create single-threaded |
| Session | `carl_setSessionLogger` | Set logger callback |
| Session | `carl_tickCallbacks` | Dispatch callbacks |
| Session | `carl_addSerializedInputSample` | Add input (bytes) |
| Session | `carl_addInputSample` | Add input (struct) |
| Session | `carl_disposeSession` | Free |
| Recognizer | `carl_createRecognizerAsync` | Create asynchronously |
| Recognizer | `carl_getCurrentScore` | Get score |
| Recognizer | `carl_setSensitivity` | Set sensitivity |
| Recognizer | `carl_getCanonicalRecordingInspector` | Get canonical inspector |
| Recognizer | `carl_disposeRecognizer` | Free |
