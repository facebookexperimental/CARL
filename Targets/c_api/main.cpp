/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/carl.h>

#ifdef CARL_PLATFORM_WINDOWS
#define C_API_EXPORT(ReturnT) __declspec(dllexport) ReturnT __cdecl
#else
#define C_API_EXPORT(ReturnT) ReturnT
#endif

extern "C"
{
    C_API_EXPORT(uint64_t) getBytes(uint64_t bytesPtr, uint8_t* destination, uint64_t size);
    C_API_EXPORT(uint64_t) startRecording();
    C_API_EXPORT(void) recordInputSample(uint64_t inProgressRecordingPtr, uint8_t* bytes, uint64_t size);
    C_API_EXPORT(uint64_t) finishRecording(uint64_t inProgressRecordingPtr);
    C_API_EXPORT(uint64_t) serializeRecording(uint64_t recordingPtr);
    C_API_EXPORT(uint64_t) deserializeRecording(uint8_t* bytes, uint64_t size);
    C_API_EXPORT(void) disposeRecording(uint64_t recordingPtr);
    C_API_EXPORT(uint64_t) getRecordingInspector(uint64_t recordingPtr);
    C_API_EXPORT(uint64_t) inspect(uint64_t recordingInspectorPtr, double timestamp);
    C_API_EXPORT(double) getStartTimestamp(uint64_t recordingInspectorPtr);
    C_API_EXPORT(double) getEndTimestamp(uint64_t recordingInspectorPtr);
    C_API_EXPORT(void) disposeRecordingInspector(uint64_t recordingInspectorPtr);
    C_API_EXPORT(uint64_t) createdDefinition(uint64_t descriptorType);
    C_API_EXPORT(void) addExample(uint64_t definitionPtr, uint64_t recordingPtr, double startTimestamp, double endTimestamp);
    C_API_EXPORT(double) getDefaultSensitivity(uint64_t definitionPtr);
    C_API_EXPORT(void) setDefaultSensitivity(uint64_t definitionPtr, double defaultSensitivity);
    C_API_EXPORT(uint64_t) serializeDefinition(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) deserializeDefinition(uint8_t* bytes, uint64_t size);
    C_API_EXPORT(uint64_t) getExamplesCount(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) getRecordingFromExampleAtIdx(uint64_t definitionPtr, uint64_t idx);
    C_API_EXPORT(double) getStartTimestampFromExampleAtIdx(uint64_t definitionPtr, uint64_t idx);
    C_API_EXPORT(double) getEndTimestampFromExampleAtIdx(uint64_t definitionPtr, uint64_t idx);
    C_API_EXPORT(void) disposeDefinition(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) createSession();
    C_API_EXPORT(void) addInputSample(uint64_t sessionPtr, uint8_t* bytes, uint64_t size);
    C_API_EXPORT(void) disposeSession(uint64_t sessionPtr);
    C_API_EXPORT(uint64_t) createRecognizer(uint64_t sessionPtr, uint64_t definitionPtr);
    C_API_EXPORT(double) getCurrentScore(uint64_t recognizerPtr);
    C_API_EXPORT(void) setSensitivity(uint64_t recognizerPtr, double sensitivity);
    C_API_EXPORT(uint64_t) getCanonicalRecordingInspector(uint64_t recognizerPtr);
    C_API_EXPORT(void) disposeRecognizer(uint64_t recognizerPtr);
}

uint64_t getBytes(uint64_t bytesPtr, uint8_t* destination, uint64_t size)
{
    auto& bytes = *reinterpret_cast<std::vector<uint8_t>*>(bytesPtr);
    if (size == 0)
    {
        return bytes.size();
    }
    else
    {
        std::memcpy(destination, bytes.data(), size);
        delete& bytes;
        return size;
    }
}

uint64_t startRecording()
{
    auto* ptr = new carl::action::InProgressRecording();
    return reinterpret_cast<uint64_t>(ptr);
}

void recordInputSample(uint64_t inProgressRecordingPtr, uint8_t* bytes, uint64_t size)
{
    auto& inProgressRecording = *reinterpret_cast<carl::action::InProgressRecording*>(inProgressRecordingPtr);
    carl::Deserialization deserialization{ bytes };
    inProgressRecording.addSample({ deserialization });
}

uint64_t finishRecording(uint64_t inProgressRecordingPtr)
{
    auto& inProgressRecording = *reinterpret_cast<carl::action::InProgressRecording*>(inProgressRecordingPtr);
    auto* recordingPtr = new carl::action::Recording(std::move(inProgressRecording));
    return reinterpret_cast<uint64_t>(recordingPtr);
}

uint64_t serializeRecording(uint64_t recordingPtr)
{
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    auto* bytesPtr = new std::vector<uint8_t>();
    carl::Serialization serialization{ *bytesPtr };
    recording.serialize(serialization);
    return reinterpret_cast<uint64_t>(bytesPtr);
}

uint64_t deserializeRecording(uint8_t* bytes, uint64_t size)
{
    carl::Deserialization deserialization{ bytes };
    auto* ptr = new carl::action::Recording(deserialization);
    return reinterpret_cast<uint64_t>(ptr);
}

void disposeRecording(uint64_t recordingPtr)
{
    delete reinterpret_cast<carl::action::Recording*>(recordingPtr);
}

uint64_t getRecordingInspector(uint64_t recordingPtr)
{
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    auto* ptr = new carl::action::RecordingInspector(recording.getInspector());
    return reinterpret_cast<uint64_t>(ptr);
}

uint64_t inspect(uint64_t recordingInspectorPtr, double timestamp)
{
    auto& inspector = *reinterpret_cast<carl::action::RecordingInspector*>(recordingInspectorPtr);
    auto& sample = inspector.inspect(timestamp);
    auto* bytesPtr = new std::vector<uint8_t>();
    carl::Serialization serialization{ *bytesPtr };
    sample.serialize(serialization);
    return reinterpret_cast<uint64_t>(bytesPtr);
}

double getStartTimestamp(uint64_t recordingInspectorPtr)
{
    auto& inspector = *reinterpret_cast<carl::action::RecordingInspector*>(recordingInspectorPtr);
    return inspector.startTimestamp();
}

double getEndTimestamp(uint64_t recordingInspectorPtr)
{
    auto& inspector = *reinterpret_cast<carl::action::RecordingInspector*>(recordingInspectorPtr);
    return inspector.endTimestamp();
}

void disposeRecordingInspector(uint64_t recordingInspectorPtr)
{
    auto* ptr = reinterpret_cast<carl::action::RecordingInspector*>(recordingInspectorPtr);
    delete ptr;
}

uint64_t createdDefinition(uint64_t descriptorType)
{
    auto* ptr = new carl::action::Definition(static_cast<carl::action::Definition::ActionType>(descriptorType));
    return reinterpret_cast<uint64_t>(ptr);
}

void addExample(
    uint64_t definitionPtr,
    uint64_t recordingPtr,
    double startTimestamp,
    double endTimestamp)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    definition.addExample({ recording, startTimestamp, endTimestamp });
}

double getDefaultSensitivity(uint64_t definitionPtr)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    return definition.DefaultSensitivity;
}

void setDefaultSensitivity(uint64_t definitionPtr, double defaultSensitivity)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    definition.DefaultSensitivity = defaultSensitivity;
}

uint64_t serializeDefinition(uint64_t definitionPtr)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto* bytesPtr = new std::vector<uint8_t>();
    carl::Serialization serialization{ *bytesPtr };
    definition.serialize(serialization);
    return reinterpret_cast<uint64_t>(bytesPtr);
}

uint64_t deserializeDefinition(uint8_t* bytes, uint64_t size)
{
    carl::Deserialization deserialization{ bytes };
    auto* ptr = new carl::action::Definition(deserialization);
    return reinterpret_cast<uint64_t>(ptr);
}

uint64_t getExamplesCount(uint64_t definitionPtr)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    return static_cast<uint64_t>(definition.getExamples().size());
}

uint64_t getRecordingFromExampleAtIdx(uint64_t definitionPtr, uint64_t idx)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto* ptr = new carl::action::Recording(definition.getExamples()[idx].getRecording());
    return reinterpret_cast<uint64_t>(ptr);
}

double getStartTimestampFromExampleAtIdx(uint64_t definitionPtr, uint64_t idx)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    return definition.getExamples()[idx].getStartTimestamp();
}

double getEndTimestampFromExampleAtIdx(uint64_t definitionPtr, uint64_t idx)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    return definition.getExamples()[idx].getEndTimestamp();
}

void disposeDefinition(uint64_t definitionPtr)
{
    delete reinterpret_cast<carl::action::Definition*>(definitionPtr);
}

uint64_t createSession()
{
    auto* ptr = new carl::Session();
    return reinterpret_cast<uint64_t>(ptr);
}

void addInputSample(uint64_t sessionPtr, uint8_t* bytes, uint64_t size)
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    carl::Deserialization deserialization{ bytes };
    session.addInput({ deserialization });
}

void disposeSession(uint64_t sessionPtr)
{
    delete reinterpret_cast<carl::Session*>(sessionPtr);
}

uint64_t createRecognizer(uint64_t sessionPtr, uint64_t definitionPtr)
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    const auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);

    auto* ptr = new carl::action::Recognizer(session, definition);
    return reinterpret_cast<uint64_t>(ptr);
}

double getCurrentScore(uint64_t recognizerPtr)
{
    auto& recognizer = *reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
    return recognizer.currentScore();
}

void setSensitivity(uint64_t recognizerPtr, double sensitivity)
{
    auto& recognizer = *reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
    recognizer.setSensitivity(sensitivity);
}

uint64_t getCanonicalRecordingInspector(uint64_t recognizerPtr)
{
    auto& recognizer = *reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
    auto* ptr = new carl::action::RecordingInspector(recognizer.getCanonicalRecordingInspector());
    return reinterpret_cast<uint64_t>(ptr);
}

void disposeRecognizer(uint64_t recognizerPtr)
{
    delete reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
}
