/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "include/carl.h"

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

    static const std::map<carl_InputSample::HAND_JOINT, carl::InputSample::Joint> openXrToCarlJointMap
    {
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_METACARPAL_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandThumb0},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_PROXIMAL_EXT, carl::InputSample::Joint::ThumbFingerBase},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_DISTAL_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandThumb2},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_TIP_EXT, carl::InputSample::Joint::ThumbFingerTip},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_PROXIMAL_EXT, carl::InputSample::Joint::IndexFingerBase},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandIndex2},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_DISTAL_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandIndex3},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_TIP_EXT, carl::InputSample::Joint::IndexFingerTip},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT, carl::InputSample::Joint::MiddleFingerBase},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandMiddle2},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_DISTAL_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandMiddle3},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_TIP_EXT, carl::InputSample::Joint::MiddleFingerTip},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_PROXIMAL_EXT, carl::InputSample::Joint::RingFingerBase},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_INTERMEDIATE_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandRing2},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_DISTAL_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandRing3},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_TIP_EXT, carl::InputSample::Joint::RingFingerTip},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_METACARPAL_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandPinky0},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_PROXIMAL_EXT, carl::InputSample::Joint::LittleFingerBase},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandPinky2},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_DISTAL_EXT, carl::InputSample::Joint::UNUSED_HandJointId_HandPinky3},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_TIP_EXT, carl::InputSample::Joint::LittleFingerTip},
    };

    carl::InputSample convert(const carl_InputSample& input)
    {
        carl::InputSample sample{};

        sample.Timestamp = input.Timestamp;

        constexpr auto cvt = [](const carl_InputSample::OptionalTransform& inT, carl::TransformT& outT) {
            outT.fromPositionOrientationScale(
                carl::VectorT
                {
                    static_cast<carl::NumberT>(inT.Position.X),
                    static_cast<carl::NumberT>(inT.Position.Y),
                    static_cast<carl::NumberT>(inT.Position.Z)
                },
                carl::QuaternionT
                {
                    static_cast<carl::NumberT>(inT.Orientation.W),
                    static_cast<carl::NumberT>(inT.Orientation.X),
                    static_cast<carl::NumberT>(inT.Orientation.Y),
                    static_cast<carl::NumberT>(inT.Orientation.Z)
                },
                carl::UNIT_SCALE);
        };

        constexpr auto cvtOptional = [cvt](const carl_InputSample::OptionalTransform& inT) {
            std::optional<carl::TransformT> outT{};
            if (inT.Valid)
            {
                outT.emplace();
                cvt(inT, *outT);
            }
            return outT;
        };

        sample.HmdPose = cvtOptional(input.HmdPose);
        sample.LeftWristPose = cvtOptional(input.LeftWristPose);
        sample.RightWristPose = cvtOptional(input.RightWristPose);
        if (input.LeftHandJointPoses[0].Valid)
        {
            sample.LeftHandJointPoses.emplace();
            for (size_t idx = 0; idx < carl_InputSample::HAND_JOINT::COUNT; ++idx)
            {
                // TODO: Contents should be replaced with the following once carl::InputSample uses the correct joints.
                // cvt(input.LeftHandJointPoses[idx], sample.LeftHandJointPoses->at(idx));
                auto found = openXrToCarlJointMap.find(static_cast<carl_InputSample::HAND_JOINT>(idx));
                if (found != openXrToCarlJointMap.end())
                {
                    cvt(input.LeftHandJointPoses[idx], sample.LeftHandJointPoses->at(static_cast<size_t>(found->second)));
                }
            }
        }
        if (input.RightHandJointPoses[0].Valid)
        {
            sample.RightHandJointPoses.emplace();
            for (size_t idx = 0; idx < carl_InputSample::HAND_JOINT::COUNT; ++idx)
            {
                // TODO: Contents should be replaced with the following once carl::InputSample uses the correct joints.
                // cvt(input.RightHandJointPoses[idx], sample.RightHandJointPoses->at(idx));
                auto found = openXrToCarlJointMap.find(static_cast<carl_InputSample::HAND_JOINT>(idx));
                if (found != openXrToCarlJointMap.end())
                {
                    cvt(input.RightHandJointPoses[idx], sample.RightHandJointPoses->at(static_cast<size_t>(found->second)));
                }
            }
        }
        
        return sample;
    }
}

uint64_t carl_getBytes(uint64_t bytesPtr, uint8_t* destination, uint64_t size)
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

uint64_t carl_startRecording(uint64_t maxSeconds)
{
    auto* ptr = new carl::action::InProgressRecording(static_cast<size_t>(maxSeconds));
    return reinterpret_cast<uint64_t>(ptr);
}

void carl_recordInputSample(uint64_t inProgressRecordingPtr, uint8_t* bytes, uint64_t size)
{
    auto& inProgressRecording = *reinterpret_cast<carl::action::InProgressRecording*>(inProgressRecordingPtr);
    carl::Deserialization deserialization{ bytes };
    inProgressRecording.addSample({ deserialization });
}

uint64_t carl_finishRecording(uint64_t inProgressRecordingPtr)
{
    auto& inProgressRecording = *reinterpret_cast<carl::action::InProgressRecording*>(inProgressRecordingPtr);
    auto* recordingPtr = new carl::action::Recording(std::move(inProgressRecording));
    return reinterpret_cast<uint64_t>(recordingPtr);
}

uint64_t carl_serializeRecording(uint64_t recordingPtr)
{
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    auto* bytesPtr = new std::vector<uint8_t>();
    carl::Serialization serialization{ *bytesPtr };
    recording.serialize(serialization);
    return reinterpret_cast<uint64_t>(bytesPtr);
}

uint64_t carl_deserializeRecording(uint8_t* bytes, uint64_t size)
{
    carl::Deserialization deserialization{ bytes };
    auto* ptr = new carl::action::Recording(deserialization);
    return reinterpret_cast<uint64_t>(ptr);
}

void carl_disposeRecording(uint64_t recordingPtr)
{
    delete reinterpret_cast<carl::action::Recording*>(recordingPtr);
}

uint64_t carl_getRecordingInspector(uint64_t recordingPtr)
{
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    auto* ptr = new carl::action::RecordingInspector(recording.getInspector());
    return reinterpret_cast<uint64_t>(ptr);
}

uint64_t carl_inspect(uint64_t recordingInspectorPtr, double timestamp)
{
    auto& inspector = *reinterpret_cast<carl::action::RecordingInspector*>(recordingInspectorPtr);
    auto& sample = inspector.inspect(timestamp);
    auto* bytesPtr = new std::vector<uint8_t>();
    carl::Serialization serialization{ *bytesPtr };
    sample.serialize(serialization);
    return reinterpret_cast<uint64_t>(bytesPtr);
}

double carl_getStartTimestamp(uint64_t recordingInspectorPtr)
{
    auto& inspector = *reinterpret_cast<carl::action::RecordingInspector*>(recordingInspectorPtr);
    return inspector.startTimestamp();
}

double carl_getEndTimestamp(uint64_t recordingInspectorPtr)
{
    auto& inspector = *reinterpret_cast<carl::action::RecordingInspector*>(recordingInspectorPtr);
    return inspector.endTimestamp();
}

void carl_disposeRecordingInspector(uint64_t recordingInspectorPtr)
{
    auto* ptr = reinterpret_cast<carl::action::RecordingInspector*>(recordingInspectorPtr);
    delete ptr;
}

uint64_t carl_createExample(uint64_t recordingPtr, double startTimestamp, double endTimestamp)
{
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    auto* examplePtr = new carl::action::Example(recording, startTimestamp, endTimestamp);
    return reinterpret_cast<uint64_t>(examplePtr);
}

uint64_t carl_createAutoTrimmedExample(uint64_t recognizerPtr, uint64_t recordingPtr)
{
    auto& recognizer = *reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    auto* examplePtr = new carl::action::Example(recognizer.createAutoTrimmedExample(recording));
    return reinterpret_cast<uint64_t>(examplePtr);
}

uint64_t carl_getRecording(uint64_t examplePtr)
{
    auto& example = *reinterpret_cast<carl::action::Example*>(examplePtr);
    auto* ptr = new carl::action::Recording(example.getRecording());
    return reinterpret_cast<uint64_t>(ptr);
}

double carl_getExampleStartTimestamp(uint64_t examplePtr)
{
    auto& example = *reinterpret_cast<carl::action::Example*>(examplePtr);
    return example.getStartTimestamp();
}

double carl_getExampleEndTimestamp(uint64_t examplePtr)
{
    auto& example = *reinterpret_cast<carl::action::Example*>(examplePtr);
    return example.getEndTimestamp();
}

void carl_disposeExample(uint64_t examplePtr)
{
    auto* ptr = reinterpret_cast<carl::action::Example*>(examplePtr);
    delete ptr;
}

uint64_t carl_createDefinition(uint64_t descriptorType)
{
    auto* ptr = new carl::action::Definition(static_cast<carl::action::Definition::ActionType>(descriptorType));
    return reinterpret_cast<uint64_t>(ptr);
}

void carl_addExample(
    uint64_t definitionPtr,
    uint64_t recordingPtr,
    double startTimestamp,
    double endTimestamp)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    definition.addExample({ recording, startTimestamp, endTimestamp });
}

void carl_addCounterexample(
    uint64_t definitionPtr,
    uint64_t recordingPtr,
    double startTimestamp,
    double endTimestamp)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto& recording = *reinterpret_cast<carl::action::Recording*>(recordingPtr);
    definition.addCounterexample({ recording, startTimestamp, endTimestamp });
}

double carl_getDefaultSensitivity(uint64_t definitionPtr)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    return definition.DefaultSensitivity;
}

void carl_setDefaultSensitivity(uint64_t definitionPtr, double defaultSensitivity)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    definition.DefaultSensitivity = defaultSensitivity;
}

uint64_t carl_serializeDefinition(uint64_t definitionPtr)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto* bytesPtr = new std::vector<uint8_t>();
    carl::Serialization serialization{ *bytesPtr };
    definition.serialize(serialization);
    return reinterpret_cast<uint64_t>(bytesPtr);
}

uint64_t carl_deserializeDefinition(uint8_t* bytes, uint64_t size)
{
    carl::Deserialization deserialization{ bytes };
    auto* ptr = new carl::action::Definition(deserialization);
    return reinterpret_cast<uint64_t>(ptr);
}

uint64_t carl_getExamplesCount(uint64_t definitionPtr)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    return static_cast<uint64_t>(definition.getExamples().size());
}

uint64_t carl_getExampleAtIdx(uint64_t definitionPtr, uint64_t idx)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto* ptr = new carl::action::Example(definition.getExamples()[idx]);
    return reinterpret_cast<uint64_t>(ptr);
}

uint64_t carl_getCounterexamplesCount(uint64_t definitionPtr)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    return static_cast<uint64_t>(definition.getCounterexamples().size());
}

uint64_t carl_getCounterexampleAtIdx(uint64_t definitionPtr, uint64_t idx)
{
    auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    auto* ptr = new carl::action::Example(definition.getCounterexamples()[idx]);
    return reinterpret_cast<uint64_t>(ptr);
}

void carl_disposeDefinition(uint64_t definitionPtr)
{
    delete reinterpret_cast<carl::action::Definition*>(definitionPtr);
}

uint64_t carl_createSession()
{
    auto* ptr = new carl::Session();
    setUpSession(*ptr);
    return reinterpret_cast<uint64_t>(ptr);
}

uint64_t carl_createSingleThreadedSession()
{
    auto* ptr = new carl::Session(true);
    setUpSession(*ptr);
    return reinterpret_cast<uint64_t>(ptr);
}

void carl_setSessionLogger(uint64_t sessionPtr, C_API_CALLBACK(void) callback(const char*))
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    session.setLogger([callback](std::string message) {
        callback(message.c_str());
    });
}

void carl_tickCallbacks(uint64_t sessionPtr)
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    session.tickCallbacks(arcana::cancellation::none());
}

void carl_addSerializedInputSample(uint64_t sessionPtr, uint8_t* bytes, uint64_t size)
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    carl::Deserialization deserialization{ bytes };
    carl::InputSample sample{ deserialization };
    session.addInput(sample);
}

void carl_addInputSample(uint64_t sessionPtr, carl_InputSample* input)
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    auto sample = convert(*input);
    session.addInput(sample);
}

void carl_disposeSession(uint64_t sessionPtr)
{
    delete reinterpret_cast<carl::Session*>(sessionPtr);
}

void carl_createRecognizerAsync(uint64_t sessionPtr, uint64_t definitionPtr, uint64_t requestId, C_API_CALLBACK(void) callback(uint64_t, uint64_t))
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    arcana::make_task(session.processingScheduler(), arcana::cancellation::none(), [&session, definitionPtr]() {
        const auto& definition = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
        return new carl::action::Recognizer(session, definition);
    }).then(session.callbackScheduler(), arcana::cancellation::none(), [callback, requestId](auto* ptr) {
        callback(requestId, reinterpret_cast<uint64_t>(ptr));
    });
}

double carl_getCurrentScore(uint64_t recognizerPtr)
{
    auto& recognizer = *reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
    return recognizer.currentScore();
}

void carl_setSensitivity(uint64_t recognizerPtr, double sensitivity)
{
    auto& recognizer = *reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
    recognizer.setSensitivity(sensitivity);
}

uint64_t carl_getCanonicalRecordingInspector(uint64_t recognizerPtr)
{
    auto& recognizer = *reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
    auto* ptr = new carl::action::RecordingInspector(recognizer.getCanonicalRecordingInspector());
    return reinterpret_cast<uint64_t>(ptr);
}

void carl_disposeRecognizer(uint64_t sessionPtr, uint64_t recognizerPtr)
{
    auto& session = *reinterpret_cast<carl::Session*>(sessionPtr);
    arcana::make_task(session.processingScheduler(), arcana::cancellation::none(), [recognizerPtr]() {
        delete reinterpret_cast<carl::action::Recognizer*>(recognizerPtr);
    });
}
