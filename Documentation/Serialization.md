# Serialization

CARL provides binary serialization for `Recording`, `Example`, and `Definition` objects. Serialized data is cross-platform and can be saved to files (conventionally using the `.carl` extension).

## C++ Serialization

### Low-Level: Serialize/Deserialization Classes

All CARL types implement `serialize()` and a deserializing constructor:

```cpp
#include <carl/Carl.h>

// Serialize a definition to bytes
std::vector<uint8_t> bytes;
carl::Serialization serialization{bytes};
definition.serialize(serialization);

// Deserialize back
carl::Deserialization deserialization{bytes.data(), bytes.size()};
carl::action::Definition restored{deserialization};
```

This pattern works for `InputSample`, `Recording`, `Example`, and `Definition`.

### High-Level: File Utilities

For most use cases, use the file serialization utilities:

```cpp
#include <carl/utilities/FileSerialization.h>

using namespace carl::utilities;
using namespace carl::action;

// Serialize to bytes
std::vector<uint8_t> bytes = Serialize(definition);

// Deserialize from bytes
auto def = TryDeserialize<Definition>(bytes);
if (def.has_value()) {
    // Use def.value()
}

// Save to file
SerializeToFile(definition, "my_gesture.carl");

// Load from file
auto loaded = TryDeserializeFromFile<Definition>("my_gesture.carl");
if (loaded.has_value()) {
    // Use loaded.value()
}
```

These utilities work with both `Definition` and `Example`.

### Legacy Format Support

If you have files from an older version of CARL:

```cpp
auto legacy = TryDeserializeLegacyFile<Definition>("old_gesture.carl");
```

## C API Serialization

The C API provides explicit serialize/deserialize functions for each type.

### Recording

```c
// Serialize
uint64_t bytesHandle = carl_serializeRecording(recordingPtr);
uint64_t size = carl_getBytes(bytesHandle, NULL, 0);
uint8_t* buffer = malloc(size);
carl_getBytes(bytesHandle, buffer, size);

// Deserialize
uint64_t recording = carl_deserializeRecording(buffer, size);
free(buffer);
```

### Definition

```c
// Serialize to bytes
uint64_t bytesHandle = carl_serializeDefinition(definitionPtr);
uint64_t size = carl_getBytes(bytesHandle, NULL, 0);
uint8_t* buffer = malloc(size);
carl_getBytes(bytesHandle, buffer, size);

// Deserialize from bytes
uint64_t def = carl_deserializeDefinition(buffer, size);
free(buffer);

// Or use file I/O directly
carl_saveDefinitionToFile(definitionPtr, "my_gesture.carl");
uint64_t loaded = carl_loadDefinitionFromFile("my_gesture.carl");
```

### Example

```c
// File I/O
carl_saveExampleToFile(examplePtr, "example_0.carl");
uint64_t loaded = carl_loadExampleFromFile("example_0.carl");
```

## Unity Serialization

### CarlDefinitionAsset (ScriptableObject)

In Unity, definitions are stored as `CarlDefinitionAsset` ScriptableObjects. These store the serialized binary data as a byte array within the Unity asset system.

**Creating a definition asset:**
1. Right-click in the Project window → **Create → CARL → Definition Asset**
2. In the Inspector, click **Import from File** and select a `.carl` file

**Programmatic access:**

```csharp
// Load a definition from an asset
CarlDefinition definition = definitionAsset.Load();
// ... use the definition ...
definition.Dispose();

// Import from file
definitionAsset.ImportFromFile("/path/to/gesture.carl");

// Set from a CarlDefinition object
definitionAsset.SetFromDefinition(definition, ActionType.RightHandGesture, "My Gesture");
```

### Programmatic Serialization in C#

```csharp
// Serialize a definition to bytes
byte[] bytes = definition.Serialize();

// Deserialize from bytes
CarlDefinition restored = CarlDefinition.Deserialize(bytes);

// File I/O
definition.SaveToFile("/path/to/gesture.carl");
CarlDefinition loaded = CarlDefinition.LoadFromFile("/path/to/gesture.carl");

// Recording serialization
byte[] recordingBytes = recording.Serialize();
CarlRecording restoredRecording = CarlRecording.Deserialize(recordingBytes);

// Example file I/O
example.SaveToFile("/path/to/example.carl");
CarlExample loadedExample = CarlExample.LoadFromFile("/path/to/example.carl");
```

Remember to call `Dispose()` on all CARL objects when done.

## File Format Notes

- CARL uses a custom binary format (not JSON, not protobuf)
- Files are cross-platform — a `.carl` file saved on Windows can be loaded on Android
- The file extension `.carl` is conventional but not enforced
- The format includes version information for forward compatibility
- Files contain all data needed to reconstruct the object, including all InputSamples in recordings
