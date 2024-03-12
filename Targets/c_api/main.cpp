/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/Carl.h>

#include <arcana/threading/task.h>

#ifdef CARL_PLATFORM_ANDROID
#include <android/log.h>
#endif

#ifdef CARL_PLATFORM_WINDOWS
#define C_API_EXPORT(ReturnT) __declspec(dllexport) ReturnT __cdecl
#define C_API_CALLBACK(ReturnT) ReturnT __cdecl
#else
#define C_API_EXPORT(ReturnT) ReturnT __cdecl
#define C_API_CALLBACK(ReturnT) ReturnT __cdecl
#endif

namespace
{
#ifdef CARL_PLATFORM_ANDROID
    void setUpSession(carl::Session& session)
    {
        session.setLogger([](std::string message) {
            __android_log_print(ANDROID_LOG_ERROR, "CARL", "%s", message.c_str());
        });
    }
#else
    void setUpSession(carl::Session& session)
    {
    }
#endif
}

extern "C"
{
    C_API_EXPORT(uint64_t) getBytes(uint64_t bytesPtr, uint8_t* destination, uint64_t size);
    C_API_EXPORT(uint64_t) startRecording(uint64_t maxSeconds);
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
    C_API_EXPORT(uint64_t) createExample(uint64_t recordingPtr, double startTimestamp, double endTimestamp);
    C_API_EXPORT(uint64_t) createAutoTrimmedExample(uint64_t recognizerPtr, uint64_t recordingPtr);
    C_API_EXPORT(uint64_t) getRecording(uint64_t examplePtr);
    C_API_EXPORT(double) getExampleStartTimestamp(uint64_t examplePtr);
    C_API_EXPORT(double) getExampleEndTimestamp(uint64_t examplePtr);
    C_API_EXPORT(void) disposeExample(uint64_t examplePtr);
    C_API_EXPORT(uint64_t) createdDefinition(uint64_t descriptorType);
    C_API_EXPORT(void) addExample(uint64_t definitionPtr, uint64_t recordingPtr, double startTimestamp, double endTimestamp);
    C_API_EXPORT(void) addCounterexample(uint64_t definitionPtr, uint64_t recordingPtr, double startTimestamp, double endTimestamp);
    C_API_EXPORT(double) getDefaultSensitivity(uint64_t definitionPtr);
    C_API_EXPORT(void) setDefaultSensitivity(uint64_t definitionPtr, double defaultSensitivity);
    C_API_EXPORT(uint64_t) serializeDefinition(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) deserializeDefinition(uint8_t* bytes, uint64_t size);
    C_API_EXPORT(uint64_t) getExamplesCount(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) getExampleAtIdx(uint64_t definitionPtr, uint64_t idx);
    C_API_EXPORT(uint64_t) getCounterexamplesCount(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) getCounterexampleAtIdx(uint64_t definitionPtr, uint64_t idx);
    C_API_EXPORT(void) disposeDefinition(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) createSession();
    C_API_EXPORT(uint64_t) createSingleThreadedSession();
    C_API_EXPORT(void) setSessionLogger(uint64_t sessionPtr, C_API_CALLBACK(void) callback(const char*));
    C_API_EXPORT(void) tickCallbacks(uint64_t sessionPtr);
    C_API_EXPORT(void) addInputSample(uint64_t sessionPtr, uint8_t* bytes, uint64_t size);
    C_API_EXPORT(void) disposeSession(uint64_t sessionPtr);
    C_API_EXPORT(void) createRecognizerAsync(uint64_t sessionPtr, uint64_t definitionPtr, uint64_t requestId, C_API_CALLBACK(void) callback(uint64_t, uint64_t));
    C_API_EXPORT(double) getCurrentScore(uint64_t recognizerPtr);
    C_API_EXPORT(void) setSensitivity(uint64_t recognizerPtr, double sensitivity);
    C_API_EXPORT(uint64_t) getCanonicalRecordingInspector(uint64_t recognizerPtr);
    C_API_EXPORT(void) disposeRecognizer(uint64_t sessionPtr, uint64_t recognizerPtr);
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

uint64_t startRecording(uint64_t maxSeconds)
{
    auto* ptr = new carl::action::InProgressRecording(static_cast<size_t>(maxSeconds));
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

uint64_t createExample(uint64_t recordingPtr, double startTimestamp, double endTimestamp)
{
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    auto* examplePtr = new carl::action::Example(recording, startTimestamp, endTimestamp);
    return reinterpret_cast<uint64_t>(examplePtr);
}

uint64_t createAutoTrimmedExample(uint64_t recognizerPtr, uint64_t recordingPtr)
{
    auto& recognizer = *reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    auto* examplePtr = new carl::action::Example(recognizer.createAutoTrimmedExample(recording));
    return reinterpret_cast<uint64_t>(examplePtr);
}

uint64_t getRecording(uint64_t examplePtr)
{
    auto& example = *reinterpret_cast<carl::action::Example*>(examplePtr);
    auto* ptr = new carl::action::Recording(example.getRecording());
    return reinterpret_cast<uint64_t>(ptr);
}

double getExampleStartTimestamp(uint64_t examplePtr)
{
    auto& example = *reinterpret_cast<carl::action::Example*>(examplePtr);
    return example.getStartTimestamp();
}

double getExampleEndTimestamp(uint64_t examplePtr)
{
    auto& example = *reinterpret_cast<carl::action::Example*>(examplePtr);
    return example.getEndTimestamp();
}

void disposeExample(uint64_t examplePtr)
{
    auto* ptr = reinterpret_cast<carl::action::Example*>(examplePtr);
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

void addCounterexample(
    uint64_t definitionPtr,
    uint64_t recordingPtr,
    double startTimestamp,
    double endTimestamp)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    definition.addCounterexample({ recording, startTimestamp, endTimestamp });
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

uint64_t getExampleAtIdx(uint64_t definitionPtr, uint64_t idx)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto* ptr = new carl::action::Example(definition.getExamples()[idx]);
    return reinterpret_cast<uint64_t>(ptr);
}

uint64_t getCounterexamplesCount(uint64_t definitionPtr)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    return static_cast<uint64_t>(definition.getCounterexamples().size());
}

uint64_t getCounterexampleAtIdx(uint64_t definitionPtr, uint64_t idx)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto* ptr = new carl::action::Example(definition.getCounterexamples()[idx]);
    return reinterpret_cast<uint64_t>(ptr);
}

void disposeDefinition(uint64_t definitionPtr)
{
    delete reinterpret_cast<carl::action::Definition*>(definitionPtr);
}

uint64_t createSession()
{
    auto* ptr = new carl::Session();
    setUpSession(*ptr);
    return reinterpret_cast<uint64_t>(ptr);
}

uint64_t createSingleThreadedSession()
{
    auto* ptr = new carl::Session(true);
    setUpSession(*ptr);
    return reinterpret_cast<uint64_t>(ptr);
}

void setSessionLogger(uint64_t sessionPtr, C_API_CALLBACK(void) callback(const char*))
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    session.setLogger([callback](std::string message) {
        callback(message.c_str());
    });
}

void tickCallbacks(uint64_t sessionPtr)
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    session.tickCallbacks(arcana::cancellation::none());
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

void createRecognizerAsync(uint64_t sessionPtr, uint64_t definitionPtr, uint64_t requestId, C_API_CALLBACK(void) callback(uint64_t, uint64_t))
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    arcana::make_task(session.processingScheduler(), arcana::cancellation::none(), [&session, definitionPtr]() {
        const auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
        return new carl::action::Recognizer(session, definition);
    }).then(session.callbackScheduler(), arcana::cancellation::none(), [callback, requestId](auto* ptr) {
        callback(requestId, reinterpret_cast<uint64_t>(ptr));
    });
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

void disposeRecognizer(uint64_t sessionPtr, uint64_t recognizerPtr)
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    arcana::make_task(session.processingScheduler(), arcana::cancellation::none(), [recognizerPtr]() {
        delete reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
    });
}
