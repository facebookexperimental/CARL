/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifdef CARL_PLATFORM_EMSCRIPTEN

#include "include/carl.h"
#include "InputSampleConversion.h"

#include <carl/Carl.h>
#include <carl/utilities/FileSerialization.h>

#include <arcana/threading/task.h>

#include <emscripten/bind.h>

namespace
{
    struct carl_InProgressRecordingWrapper
    {
        void addInput(const carl_InputSample& capiSample)
        {
            InProgressRecording.addSample(carl::capi::convert(capiSample));
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
            return carl::capi::convert(Inspector.inspect(timestamp));
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

        size_t getActionType() const
        {
            return static_cast<size_t>(Definition.getDescriptorType());
        }

        void addExample(const carl_ExampleWrapper& example)
        {
            Definition.addExample(example.Example);
        }

        size_t getExamplesCount() const
        {
            return Definition.getExamples().size();
        }

        carl_ExampleWrapper getExample(size_t idx) const
        {
            auto example{Definition.getExamples()[idx]};
            return{std::move(example)};
        }

        void addCounterexample(const carl_ExampleWrapper& counterexample)
        {
            Definition.addCounterexample(counterexample.Example);
        }

        size_t getCounterexamplesCount() const
        {
            return Definition.getCounterexamples().size();
        }

        carl_ExampleWrapper getCounterexample(size_t idx) const
        {
            auto counterexample{Definition.getCounterexamples()[idx]};
            return{std::move(counterexample)};
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
            auto sample = carl::capi::convert(capiSample);
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
}

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
        .function("getActionType", &carl_DefinitionWrapper::getActionType)
        .function("addExample", &carl_DefinitionWrapper::addExample)
        .function("getExamplesCount", &carl_DefinitionWrapper::getExamplesCount)
        .function("getExample", &carl_DefinitionWrapper::getExample)
        .function("addCounterexample", &carl_DefinitionWrapper::addCounterexample)
        .function("getCounterexamplesCount", &carl_DefinitionWrapper::getCounterexamplesCount)
        .function("getCounterexample", &carl_DefinitionWrapper::getCounterexample)
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
