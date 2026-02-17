/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "include/carl.h"

#include <carl/Carl.h>
#include <carl/utilities/FileSerialization.h>

#include <arcana/threading/task.h>

#ifdef CARL_PLATFORM_ANDROID
#include <android/log.h>
#endif

#ifdef CARL_PLATFORM_EMSCRIPTEN
#include <emscripten/bind.h>
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
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_METACARPAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_THUMB_METACARPAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_PROXIMAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_THUMB_PROXIMAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_DISTAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_THUMB_DISTAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_TIP_EXT, carl::InputSample::Joint::XR_HAND_JOINT_THUMB_TIP_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_PROXIMAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_INDEX_PROXIMAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT, carl::InputSample::Joint::XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_DISTAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_INDEX_DISTAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_TIP_EXT, carl::InputSample::Joint::XR_HAND_JOINT_INDEX_TIP_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT, carl::InputSample::Joint::XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_DISTAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_MIDDLE_DISTAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_TIP_EXT, carl::InputSample::Joint::XR_HAND_JOINT_MIDDLE_TIP_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_PROXIMAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_RING_PROXIMAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_INTERMEDIATE_EXT, carl::InputSample::Joint::XR_HAND_JOINT_RING_INTERMEDIATE_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_DISTAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_RING_DISTAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_TIP_EXT, carl::InputSample::Joint::XR_HAND_JOINT_RING_TIP_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_METACARPAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_LITTLE_METACARPAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_PROXIMAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_LITTLE_PROXIMAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT, carl::InputSample::Joint::XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_DISTAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_LITTLE_DISTAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_TIP_EXT, carl::InputSample::Joint::XR_HAND_JOINT_LITTLE_TIP_EXT},
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
        if (input.LeftControllerState.Valid)
        {
            sample.LeftControllerInput = 
                std::array<carl::NumberT, static_cast<size_t>(carl::InputSample::ControllerInput::COUNT)>{
                static_cast<carl::NumberT>(input.LeftControllerState.PrimaryClick),
                static_cast<carl::NumberT>(input.LeftControllerState.SecondaryClick),
                static_cast<carl::NumberT>(input.LeftControllerState.ThumbstickX),
                static_cast<carl::NumberT>(input.LeftControllerState.ThumbstickY),
                static_cast<carl::NumberT>(input.LeftControllerState.ThumbstickClick),
                static_cast<carl::NumberT>(input.LeftControllerState.SqueezeValue),
                static_cast<carl::NumberT>(input.LeftControllerState.TriggerValue),
            };
        }
        if (input.RightControllerState.Valid)
        {
            sample.RightControllerInput =
                std::array<carl::NumberT, static_cast<size_t>(carl::InputSample::ControllerInput::COUNT)>{
                static_cast<carl::NumberT>(input.RightControllerState.PrimaryClick),
                static_cast<carl::NumberT>(input.RightControllerState.SecondaryClick),
                static_cast<carl::NumberT>(input.RightControllerState.ThumbstickX),
                static_cast<carl::NumberT>(input.RightControllerState.ThumbstickY),
                static_cast<carl::NumberT>(input.RightControllerState.ThumbstickClick),
                static_cast<carl::NumberT>(input.RightControllerState.SqueezeValue),
                static_cast<carl::NumberT>(input.RightControllerState.TriggerValue),
            };
        }
        
        return sample;
    }

    carl_InputSample convert(const carl::InputSample& input)
    {
        carl_InputSample sample{};

        sample.Timestamp = input.Timestamp;

        constexpr auto cvt = [](const carl::TransformT& inT, carl_InputSample::OptionalTransform& outT) {
            carl::VectorT position{ inT.translation() };
            carl::QuaternionT orientation{ inT.rotation() };
            outT.Position.X = position.x();
            outT.Position.Y = position.y();
            outT.Position.Z = position.z();
            outT.Orientation.W = orientation.w();
            outT.Orientation.X = orientation.x();
            outT.Orientation.Y = orientation.y();
            outT.Orientation.Z = orientation.z();
        };

        constexpr auto cvtOptional = [cvt](const std::optional<carl::TransformT>& inT) {
            carl_InputSample::OptionalTransform outT{};
            if (inT.has_value())
            {
                outT.Valid = true;
                cvt(*inT, outT);
            }
            return outT;
        };

        sample.HmdPose = cvtOptional(input.HmdPose);
        sample.LeftWristPose = cvtOptional(input.LeftWristPose);
        sample.RightWristPose = cvtOptional(input.RightWristPose);
        if (input.LeftHandJointPoses.has_value())
        {
            for (size_t idx = 0; idx < carl_InputSample::HAND_JOINT::COUNT; ++idx)
            {
                // TODO: Contents should be replaced with the following once carl::InputSample uses the correct joints.
                // cvt(input.LeftHandJointPoses[idx], sample.LeftHandJointPoses->at(idx));
                auto found = openXrToCarlJointMap.find(static_cast<carl_InputSample::HAND_JOINT>(idx));
                if (found != openXrToCarlJointMap.end())
                {
                    cvt(input.LeftHandJointPoses->at(idx), sample.LeftHandJointPoses[static_cast<size_t>(found->second)]);
                }
            }
        }
        if (input.RightHandJointPoses.has_value())
        {
            for (size_t idx = 0; idx < carl_InputSample::HAND_JOINT::COUNT; ++idx)
            {
                // TODO: Contents should be replaced with the following once carl::InputSample uses the correct joints.
                // cvt(input.RightHandJointPoses[idx], sample.RightHandJointPoses->at(idx));
                auto found = openXrToCarlJointMap.find(static_cast<carl_InputSample::HAND_JOINT>(idx));
                if (found != openXrToCarlJointMap.end())
                {
                    cvt(input.RightHandJointPoses->at(idx), sample.RightHandJointPoses[static_cast<size_t>(found->second)]);
                }
            }
        }
        if (input.LeftControllerInput.has_value())
        {
            sample.LeftControllerState.Valid = true;
            sample.LeftControllerState.PrimaryClick = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_PRIMARY_CLICK));
            sample.LeftControllerState.SecondaryClick = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_SECONDARY_CLICK));
            sample.LeftControllerState.ThumbstickX = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_X));
            sample.LeftControllerState.ThumbstickY = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_Y));
            sample.LeftControllerState.ThumbstickClick = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_CLICK));
            sample.LeftControllerState.SqueezeValue = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_SQUEEZE_VALUE));
            sample.LeftControllerState.TriggerValue = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_TRIGGER_VALUE));
        }
        if (input.RightControllerInput.has_value())
        {
            sample.RightControllerState.Valid = true;
            sample.RightControllerState.PrimaryClick = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_PRIMARY_CLICK));
            sample.RightControllerState.SecondaryClick = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_SECONDARY_CLICK));
            sample.RightControllerState.ThumbstickX = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_X));
            sample.RightControllerState.ThumbstickY = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_Y));
            sample.RightControllerState.ThumbstickClick = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_CLICK));
            sample.RightControllerState.SqueezeValue = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_SQUEEZE_VALUE));
            sample.RightControllerState.TriggerValue = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_TRIGGER_VALUE));
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
        std::memcpy(destination, bytes.data(), static_cast<size_t>(size));
        delete& bytes;
        return size;
    }
}

uint64_t carl_startRecording(uint64_t maxSeconds)
{
    auto* ptr = new carl::action::InProgressRecording(static_cast<size_t>(maxSeconds));
    return reinterpret_cast<uint64_t>(ptr);
}

void carl_recordObjectInputSample(uint64_t inProgressRecordingPtr, carl_InputSample* sample)
{
    auto& inProgressRecording = *reinterpret_cast<carl::action::InProgressRecording*>(inProgressRecordingPtr);
    inProgressRecording.addSample(convert(*sample));
}

void carl_recordInputSample(uint64_t inProgressRecordingPtr, uint8_t* bytes, uint64_t size)
{
    auto& inProgressRecording = *reinterpret_cast<carl::action::InProgressRecording*>(inProgressRecordingPtr);
    carl::Deserialization deserialization{ bytes, static_cast<size_t>(size) };
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
    carl::Deserialization deserialization{ bytes, static_cast<size_t>(size) };
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
    auto* ptr = new carl::action::Definition(static_cast<carl::action::ActionType>(descriptorType));
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
    carl::Deserialization deserialization{ bytes, static_cast<size_t>(size) };
    auto* ptr = new carl::action::Definition(deserialization);
    return reinterpret_cast<uint64_t>(ptr);
}

uint64_t carl_loadExampleFromFile(const char* path)
{
    auto loaded = carl::utilities::TryDeserializeFromFile<carl::action::Example>(path);
    if (!loaded.has_value())
    {
        return 0;
    }
    auto* ptr = new carl::action::Example(std::move(loaded.value()));
    return reinterpret_cast<uint64_t>(ptr);
}

void carl_saveExampleToFile(uint64_t examplePtr, const char* path)
{
    const auto& example = *reinterpret_cast<carl::action::Example*>(examplePtr);
    carl::utilities::SerializeToFile<carl::action::Example>(example, path);
}

uint64_t carl_loadDefinitionFromFile(const char* path)
{
    auto loaded = carl::utilities::TryDeserializeFromFile<carl::action::Definition>(path);
    if (!loaded.has_value())
    {
        return 0;
    }
    auto* ptr = new carl::action::Definition(std::move(loaded.value()));
    return reinterpret_cast<uint64_t>(ptr);
}

void carl_saveDefinitionToFile(uint64_t definitionPtr, const char* path)
{
    const auto& example = *reinterpret_cast<carl::action::Definition*>(definitionPtr);
    carl::utilities::SerializeToFile<carl::action::Definition>(example, path);
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
    carl::Deserialization deserialization{ bytes, static_cast<size_t>(size) };
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

#ifdef CARL_PLATFORM_EMSCRIPTEN
struct carl_InProgressRecordingWrapper
{
    void addInput(const carl_InputSample& capiSample)
    {
        InProgressRecording.addSample(convert(capiSample));
    }

    carl::action::InProgressRecording InProgressRecording{};
};

struct carl_RecordingWrapper
{
    carl_RecordingWrapper(carl::action::Recording recording)
        : Recording{ std::move(recording) }
    {
    }

    carl_RecordingWrapper(carl_InProgressRecordingWrapper& wrapper)
        : Recording{ std::move(wrapper.InProgressRecording) }
    {
    }

    carl::action::Recording Recording;
};

struct carl_RecordingInspectorWrapper
{
    carl_RecordingInspectorWrapper(const carl_RecordingWrapper& recording)
        : Inspector{ recording.Recording.getInspector() }
    {
    }

    carl_RecordingInspectorWrapper(carl::action::RecordingInspector inspector)
        : Inspector{ std::move(inspector) }
    {
    }

    double startTimestamp() const
    {
        return Inspector.startTimestamp();
    }

    double endTimestamp() const
    {
        return Inspector.endTimestamp();
    }

    carl_InputSample inspect(double timestamp)
    {
        return convert(Inspector.inspect(timestamp));
    }

    carl::action::RecordingInspector Inspector;
};

struct carl_ExampleWrapper
{
    carl_ExampleWrapper(const carl_RecordingWrapper& recording, double startTimestamp, double endTimestamp)
        : Example{ recording.Recording, startTimestamp, endTimestamp }
        , RecordingWrapper{ recording }
    {
    }
    
    carl_ExampleWrapper(carl::action::Example&& example)
        : Example{ std::forward<carl::action::Example>(example) }
        , RecordingWrapper{ Example.getRecording() }
    {
    }

    const carl_RecordingWrapper& getRecording() const
    {
        return RecordingWrapper;
    }

    carl_RecordingInspectorWrapper getRecordingInspector() const
    {
        return Example.getRecording().getInspector();
    }

    double getStartTimestamp() const
    {
        return Example.getStartTimestamp();
    }

    void setStartTimestamp(double timestamp)
    {
        Example.setStartTimestamp(timestamp);
    }

    double getEndTimestamp() const
    {
        return Example.getEndTimestamp();
    }

    void setEndTimestamp(double timestamp)
    {
        Example.setEndTimestamp(timestamp);
    }

    std::vector<uint8_t> serialize() const
    {
        return carl::utilities::Serialize(Example);
    }

    static std::optional<carl_ExampleWrapper> tryDeserialize(const std::vector<uint8_t>& bytes)
    {
        auto example = carl::utilities::TryDeserialize<carl::action::Example>(bytes);
        if (example.has_value())
        {
            return carl_ExampleWrapper{ std::move(example.value()) };
        }
        else
        {
            return{};
        }
    }

    carl::action::Example Example;
    carl_RecordingWrapper RecordingWrapper;
};

// Emscripten does not allow for 64-bit enums, so we have to wrap this.
enum class carl_ActionTypeWrapper : uint32_t
{
    LeftHandPose = static_cast<uint32_t>(carl::action::ActionType::LeftHandPose),
    LeftHandGesture = static_cast<uint32_t>(carl::action::ActionType::LeftHandGesture),
    RightHandPose = static_cast<uint32_t>(carl::action::ActionType::RightHandPose),
    RightHandGesture = static_cast<uint32_t>(carl::action::ActionType::RightHandGesture),
    TwoHandGesture = static_cast<uint32_t>(carl::action::ActionType::TwoHandGesture),
    LeftControllerGesture = static_cast<uint32_t>(carl::action::ActionType::LeftControllerGesture),
    RightControllerGesture = static_cast<uint32_t>(carl::action::ActionType::RightControllerGesture),
    TwoControllerGesture = static_cast<uint32_t>(carl::action::ActionType::TwoControllerGesture),
    LeftWristTrajectory = static_cast<uint32_t>(carl::action::ActionType::LeftWristTrajectory),
    RightWristTrajectory = static_cast<uint32_t>(carl::action::ActionType::RightWristTrajectory),
    LeftHandShape = static_cast<uint32_t>(carl::action::ActionType::LeftHandShape),
    RightHandShape = static_cast<uint32_t>(carl::action::ActionType::RightHandShape),
};

struct carl_DefinitionWrapper
{
    carl_DefinitionWrapper(carl_ActionTypeWrapper type)
        : Definition{ static_cast<carl::action::ActionType>(type) }
    {
    }

    carl_DefinitionWrapper(carl::action::Definition&& definition)
        : Definition{ std::forward<carl::action::Definition>(definition) }
    {
    }

    void addExample(const carl_ExampleWrapper& example)
    {
        Definition.addExample(example.Example);
    }

    void addCounterexample(const carl_ExampleWrapper& counterexample)
    {
        Definition.addCounterexample(counterexample.Example);
    }

    auto getDefaultSensitivity() const
    {
        return Definition.DefaultSensitivity;
    }

    void setDefaultSensitivity(double sensitivity)
    {
        Definition.DefaultSensitivity = sensitivity;
    }

    std::vector<uint8_t> serialize() const
    {
        return carl::utilities::Serialize(Definition);
    }

    static std::optional<carl_DefinitionWrapper> tryDeserialize(const std::vector<uint8_t>& bytes)
    {
        auto definition = carl::utilities::TryDeserialize<carl::action::Definition>(bytes);
        if (definition.has_value())
        {
            return carl_DefinitionWrapper{ std::move(definition.value()) };
        }
        else
        {
            return{};
        }
    }

    carl::action::Definition Definition;
};

struct carl_SessionWrapper
{
    void addInput(const carl_InputSample& capiSample)
    {
        auto sample = convert(capiSample);
        Session.addInput(sample);
    }

    void tickCallbacks()
    {
        Session.tickCallbacks(arcana::cancellation::none());
    }

    carl::Session Session{ true };
};

struct carl_RecognizerWrapper
{
    carl_RecognizerWrapper(carl_SessionWrapper& sessionWrapper, const carl_DefinitionWrapper& definitionWrapper)
        : Recognizer{}
    {
        auto& session = sessionWrapper.Session;
        auto& definition = definitionWrapper.Definition;
        Recognizer = std::make_unique<carl::action::Recognizer>(session, definition);
        // TODO: This will need to be made asynchronous if we want to support multithreading.
        /*arcana::make_task(session.processingScheduler(), arcana::cancellation::none(), [&session, &definition]() {
            return std::make_unique<carl::action::Recognizer>(session, definition);
        }).then(session.callbackScheduler(), arcana::cancellation::none(), [this](std::unique_ptr<carl::action::Recognizer>& ptr) mutable {
            Recognizer.swap(ptr);
        });*/
    }

    double currentScore()
    {
        if (Recognizer != nullptr)
        {
            return Recognizer->currentScore();
        }
        else
        {
            return 0;
        }
    }

    carl_RecordingInspectorWrapper getCanonicalRecordingInspector() const
    {
        auto inspector = Recognizer->getCanonicalRecordingInspector();
        return{ std::move(inspector) };
    }

    double getSensitivity() const
    {
        return Recognizer->getSensitivity();
    }

    void setSensitivity(double sensitivity)
    {
        if (Recognizer != nullptr)
        {
            Recognizer->setSensitivity(sensitivity);
        }
    }

    std::unique_ptr<carl::action::Recognizer> Recognizer{};
};

carl_InputSample carl_createInputSample() { return{}; }

EMSCRIPTEN_BINDINGS(carl_bindings) {
    emscripten::register_vector<uint8_t>("SerializedBytes");

    emscripten::enum_<carl_InputSample::HAND_JOINT>("HAND_JOINT")
        .value("XR_HAND_JOINT_PALM_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_PALM_EXT)
        .value("XR_HAND_JOINT_WRIST_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_WRIST_EXT)
        .value("XR_HAND_JOINT_THUMB_METACARPAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_METACARPAL_EXT)
        .value("XR_HAND_JOINT_THUMB_PROXIMAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_PROXIMAL_EXT)
        .value("XR_HAND_JOINT_THUMB_DISTAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_DISTAL_EXT)
        .value("XR_HAND_JOINT_THUMB_TIP_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_TIP_EXT)
        .value("XR_HAND_JOINT_INDEX_METACARPAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_METACARPAL_EXT)
        .value("XR_HAND_JOINT_INDEX_PROXIMAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_PROXIMAL_EXT)
        .value("XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT)
        .value("XR_HAND_JOINT_INDEX_DISTAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_DISTAL_EXT)
        .value("XR_HAND_JOINT_INDEX_TIP_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_TIP_EXT)
        .value("XR_HAND_JOINT_MIDDLE_METACARPAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_METACARPAL_EXT)
        .value("XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT)
        .value("XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT)
        .value("XR_HAND_JOINT_MIDDLE_DISTAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_DISTAL_EXT)
        .value("XR_HAND_JOINT_MIDDLE_TIP_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_TIP_EXT)
        .value("XR_HAND_JOINT_RING_METACARPAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_METACARPAL_EXT)
        .value("XR_HAND_JOINT_RING_PROXIMAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_PROXIMAL_EXT)
        .value("XR_HAND_JOINT_RING_INTERMEDIATE_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_INTERMEDIATE_EXT)
        .value("XR_HAND_JOINT_RING_DISTAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_DISTAL_EXT)
        .value("XR_HAND_JOINT_RING_TIP_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_TIP_EXT)
        .value("XR_HAND_JOINT_LITTLE_METACARPAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_METACARPAL_EXT)
        .value("XR_HAND_JOINT_LITTLE_PROXIMAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_PROXIMAL_EXT)
        .value("XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT)
        .value("XR_HAND_JOINT_LITTLE_DISTAL_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_DISTAL_EXT)
        .value("XR_HAND_JOINT_LITTLE_TIP_EXT", carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_TIP_EXT)
        .value("COUNT", carl_InputSample::HAND_JOINT::COUNT);

    emscripten::value_object<decltype(carl_InputSample::OptionalTransform::Position)>("Position")
        .field("x", &decltype(carl_InputSample::OptionalTransform::Position)::X)
        .field("y", &decltype(carl_InputSample::OptionalTransform::Position)::Y)
        .field("z", &decltype(carl_InputSample::OptionalTransform::Position)::Z);

    emscripten::value_object<decltype(carl_InputSample::OptionalTransform::Orientation)>("Orientation")
        .field("w", &decltype(carl_InputSample::OptionalTransform::Orientation)::W)
        .field("x", &decltype(carl_InputSample::OptionalTransform::Orientation)::X)
        .field("y", &decltype(carl_InputSample::OptionalTransform::Orientation)::Y)
        .field("z", &decltype(carl_InputSample::OptionalTransform::Orientation)::Z);

    emscripten::value_object<carl_InputSample::OptionalTransform>("OptionalTransform")
        .field("valid", &carl_InputSample::OptionalTransform::Valid)
        .field("position", &carl_InputSample::OptionalTransform::Position)
        .field("orientation", &carl_InputSample::OptionalTransform::Orientation);

    emscripten::value_array<std::array<carl_InputSample::OptionalTransform, 26>>("OptionalTransformArray")
        .element(emscripten::index<0>())
        .element(emscripten::index<1>())
        .element(emscripten::index<2>())
        .element(emscripten::index<3>())
        .element(emscripten::index<4>())
        .element(emscripten::index<5>())
        .element(emscripten::index<6>())
        .element(emscripten::index<7>())
        .element(emscripten::index<8>())
        .element(emscripten::index<9>())
        .element(emscripten::index<10>())
        .element(emscripten::index<11>())
        .element(emscripten::index<12>())
        .element(emscripten::index<13>())
        .element(emscripten::index<14>())
        .element(emscripten::index<15>())
        .element(emscripten::index<16>())
        .element(emscripten::index<17>())
        .element(emscripten::index<18>())
        .element(emscripten::index<19>())
        .element(emscripten::index<20>())
        .element(emscripten::index<21>())
        .element(emscripten::index<22>())
        .element(emscripten::index<23>())
        .element(emscripten::index<24>())
        .element(emscripten::index<25>());

    emscripten::value_object<carl_InputSample::OptionalControllerState>("OptionalControllerState")
        .field("valid", &carl_InputSample::OptionalControllerState::Valid)
        .field("primaryClick", &carl_InputSample::OptionalControllerState::PrimaryClick)
        .field("secondaryClick", &carl_InputSample::OptionalControllerState::SecondaryClick)
        .field("thumbstickX", &carl_InputSample::OptionalControllerState::ThumbstickX)
        .field("thumbstickY", &carl_InputSample::OptionalControllerState::ThumbstickY)
        .field("thumbstickClick", &carl_InputSample::OptionalControllerState::ThumbstickClick)
        .field("squeezeValue", &carl_InputSample::OptionalControllerState::SqueezeValue)
        .field("triggerValue", &carl_InputSample::OptionalControllerState::TriggerValue);

    emscripten::value_object<carl_InputSample>("InputSample")
        .field("timestamp", &carl_InputSample::Timestamp)
        .field("hmdPose", &carl_InputSample::HmdPose)
        .field("leftWristPose", &carl_InputSample::LeftWristPose)
        .field("rightWristPose", &carl_InputSample::RightWristPose)
        .field("leftHandJointPoses", &carl_InputSample::LeftHandJointPoses)
        .field("rightHandJointPoses", &carl_InputSample::RightHandJointPoses)
        .field("leftControllerState", &carl_InputSample::LeftControllerState)
        .field("rightControllerState", &carl_InputSample::RightControllerState);

    emscripten::function("createInputSample", &carl_createInputSample, emscripten::return_value_policy::take_ownership());

    emscripten::class_<carl_InProgressRecordingWrapper>("InProgressRecording")
        .constructor<>()
        .function("addInput", &carl_InProgressRecordingWrapper::addInput);

    emscripten::class_<carl_RecordingWrapper>("Recording")
        .constructor<carl_InProgressRecordingWrapper&>();

    emscripten::class_<carl_RecordingInspectorWrapper>("RecordingInspector")
        .constructor<carl_RecordingWrapper&>()
        .function("startTimestamp", &carl_RecordingInspectorWrapper::startTimestamp)
        .function("endTimestamp", &carl_RecordingInspectorWrapper::endTimestamp)
        .function("inspect", &carl_RecordingInspectorWrapper::inspect);

    emscripten::register_optional<carl_ExampleWrapper>();
    emscripten::class_<carl_ExampleWrapper>("Example")
        .constructor<carl_RecordingWrapper&, double, double>()
        .function("getRecording", &carl_ExampleWrapper::getRecording)
        .function("getRecordingInspector", &carl_ExampleWrapper::getRecordingInspector)
        .function("getStartTimestamp", &carl_ExampleWrapper::getStartTimestamp)
        .function("setStartTimestamp", &carl_ExampleWrapper::setStartTimestamp)
        .function("getEndTimestamp", &carl_ExampleWrapper::getEndTimestamp)
        .function("setEndTimestamp", &carl_ExampleWrapper::setEndTimestamp)
        .function("serialize", &carl_ExampleWrapper::serialize)
        .class_function("tryDeserialize", &carl_ExampleWrapper::tryDeserialize);

    emscripten::enum_<carl_ActionTypeWrapper>("ACTION_TYPE")
        .value("LeftHandPose", carl_ActionTypeWrapper::LeftHandPose)
        .value("LeftHandGesture", carl_ActionTypeWrapper::LeftHandGesture)
        .value("RightHandPose", carl_ActionTypeWrapper::RightHandPose)
        .value("RightHandGesture", carl_ActionTypeWrapper::RightHandGesture)
        .value("TwoHandGesture", carl_ActionTypeWrapper::TwoHandGesture)
        .value("LeftControllerGesture", carl_ActionTypeWrapper::LeftControllerGesture)
        .value("RightControllerGesture", carl_ActionTypeWrapper::RightControllerGesture)
        .value("TwoControllerGesture", carl_ActionTypeWrapper::TwoControllerGesture)
        .value("LeftWristTrajectory", carl_ActionTypeWrapper::LeftWristTrajectory)
        .value("RightWristTrajectory", carl_ActionTypeWrapper::RightWristTrajectory)
        .value("LeftHandShape", carl_ActionTypeWrapper::LeftHandShape)
        .value("RightHandShape", carl_ActionTypeWrapper::RightHandShape);

    emscripten::register_optional<carl_DefinitionWrapper>();
    emscripten::class_<carl_DefinitionWrapper>("Definition")
        .constructor<carl_ActionTypeWrapper>()
        .function("addExample", &carl_DefinitionWrapper::addExample)
        .function("addCounterexample", &carl_DefinitionWrapper::addCounterexample)
        .function("getDefaultSensitivity", &carl_DefinitionWrapper::getDefaultSensitivity)
        .function("setDefaultSensitivity", &carl_DefinitionWrapper::setDefaultSensitivity)
        .function("serialize", &carl_DefinitionWrapper::serialize)
        .class_function("tryDeserialize", &carl_DefinitionWrapper::tryDeserialize);

    emscripten::class_<carl_SessionWrapper>("Session")
        .constructor<>()
        .function("addInput", &carl_SessionWrapper::addInput)
        .function("tickCallbacks", &carl_SessionWrapper::tickCallbacks);

    emscripten::class_<carl_RecognizerWrapper>("Recognizer")
        .constructor<carl_SessionWrapper&, carl_DefinitionWrapper&>()
        .function("currentScore", &carl_RecognizerWrapper::currentScore)
        .function("getCanonicalRecordingInspector", &carl_RecognizerWrapper::getCanonicalRecordingInspector)
        .function("getSensitivity", &carl_RecognizerWrapper::getSensitivity)
        .function("setSensitivity", &carl_RecognizerWrapper::setSensitivity);
}
#endif
