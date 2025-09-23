#pragma once

#include <cstdint>

#ifdef CARL_PLATFORM_WINDOWS
#define C_API_EXPORT(ReturnT) __declspec(dllexport) ReturnT __cdecl
#define C_API_CALLBACK(ReturnT) ReturnT __cdecl
#else
#define C_API_EXPORT(ReturnT) ReturnT __cdecl
#define C_API_CALLBACK(ReturnT) ReturnT __cdecl
#endif

extern "C"
{
    struct carl_InputSample
    {
        enum HAND_JOINT
        {
            XR_HAND_JOINT_PALM_EXT = 0,
            XR_HAND_JOINT_WRIST_EXT = 1,
            XR_HAND_JOINT_THUMB_METACARPAL_EXT = 2,
            XR_HAND_JOINT_THUMB_PROXIMAL_EXT = 3,
            XR_HAND_JOINT_THUMB_DISTAL_EXT = 4,
            XR_HAND_JOINT_THUMB_TIP_EXT = 5,
            XR_HAND_JOINT_INDEX_METACARPAL_EXT = 6,
            XR_HAND_JOINT_INDEX_PROXIMAL_EXT = 7,
            XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT = 8,
            XR_HAND_JOINT_INDEX_DISTAL_EXT = 9,
            XR_HAND_JOINT_INDEX_TIP_EXT = 10,
            XR_HAND_JOINT_MIDDLE_METACARPAL_EXT = 11,
            XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT = 12,
            XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT = 13,
            XR_HAND_JOINT_MIDDLE_DISTAL_EXT = 14,
            XR_HAND_JOINT_MIDDLE_TIP_EXT = 15,
            XR_HAND_JOINT_RING_METACARPAL_EXT = 16,
            XR_HAND_JOINT_RING_PROXIMAL_EXT = 17,
            XR_HAND_JOINT_RING_INTERMEDIATE_EXT = 18,
            XR_HAND_JOINT_RING_DISTAL_EXT = 19,
            XR_HAND_JOINT_RING_TIP_EXT = 20,
            XR_HAND_JOINT_LITTLE_METACARPAL_EXT = 21,
            XR_HAND_JOINT_LITTLE_PROXIMAL_EXT = 22,
            XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT = 23,
            XR_HAND_JOINT_LITTLE_DISTAL_EXT = 24,
            XR_HAND_JOINT_LITTLE_TIP_EXT = 25,
            COUNT = 26,
        };

        struct OptionalTransform
        {
            bool Valid{};

            struct
            {
                double X{}, Y{}, Z{};
            } Position{};

            struct
            {
                double W{}, X{}, Y{}, Z{};
            } Orientation{};
        };

        double Timestamp{};

        OptionalTransform HmdPose{};
        OptionalTransform LeftWristPose{};
        OptionalTransform RightWristPose{};
        OptionalTransform LeftHandJointPoses[HAND_JOINT::COUNT];
        OptionalTransform RightHandJointPoses[HAND_JOINT::COUNT];
    };

    C_API_EXPORT(uint64_t) carl_getBytes(uint64_t bytesPtr, uint8_t* destination, uint64_t size);
    C_API_EXPORT(uint64_t) carl_startRecording(uint64_t maxSeconds);
    C_API_EXPORT(void) carl_recordInputSample(uint64_t inProgressRecordingPtr, uint8_t* bytes, uint64_t size);
    C_API_EXPORT(uint64_t) carl_finishRecording(uint64_t inProgressRecordingPtr);
    C_API_EXPORT(uint64_t) carl_serializeRecording(uint64_t recordingPtr);
    C_API_EXPORT(uint64_t) carl_deserializeRecording(uint8_t* bytes, uint64_t size);
    C_API_EXPORT(void) carl_disposeRecording(uint64_t recordingPtr);
    C_API_EXPORT(uint64_t) carl_getRecordingInspector(uint64_t recordingPtr);
    C_API_EXPORT(uint64_t) carl_inspect(uint64_t recordingInspectorPtr, double timestamp);
    C_API_EXPORT(double) carl_getStartTimestamp(uint64_t recordingInspectorPtr);
    C_API_EXPORT(double) carl_getEndTimestamp(uint64_t recordingInspectorPtr);
    C_API_EXPORT(void) carl_disposeRecordingInspector(uint64_t recordingInspectorPtr);
    C_API_EXPORT(uint64_t) carl_createExample(uint64_t recordingPtr, double startTimestamp, double endTimestamp);
    C_API_EXPORT(uint64_t) carl_createAutoTrimmedExample(uint64_t recognizerPtr, uint64_t recordingPtr);
    C_API_EXPORT(uint64_t) carl_getRecording(uint64_t examplePtr);
    C_API_EXPORT(double) carl_getExampleStartTimestamp(uint64_t examplePtr);
    C_API_EXPORT(double) carl_getExampleEndTimestamp(uint64_t examplePtr);
    C_API_EXPORT(void) carl_disposeExample(uint64_t examplePtr);
    C_API_EXPORT(uint64_t) carl_createDefinition(uint64_t descriptorType);
    C_API_EXPORT(void) carl_addExample(uint64_t definitionPtr, uint64_t recordingPtr, double startTimestamp, double endTimestamp);
    C_API_EXPORT(void) carl_addCounterexample(uint64_t definitionPtr, uint64_t recordingPtr, double startTimestamp, double endTimestamp);
    C_API_EXPORT(double) carl_getDefaultSensitivity(uint64_t definitionPtr);
    C_API_EXPORT(void) carl_setDefaultSensitivity(uint64_t definitionPtr, double defaultSensitivity);
    C_API_EXPORT(uint64_t) carl_serializeDefinition(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) carl_deserializeDefinition(uint8_t* bytes, uint64_t size);
    C_API_EXPORT(uint64_t) carl_getExamplesCount(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) carl_getExampleAtIdx(uint64_t definitionPtr, uint64_t idx);
    C_API_EXPORT(uint64_t) carl_getCounterexamplesCount(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) carl_getCounterexampleAtIdx(uint64_t definitionPtr, uint64_t idx);
    C_API_EXPORT(void) carl_disposeDefinition(uint64_t definitionPtr);
    C_API_EXPORT(uint64_t) carl_createSession();
    C_API_EXPORT(uint64_t) carl_createSingleThreadedSession();
    C_API_EXPORT(void) carl_setSessionLogger(uint64_t sessionPtr, C_API_CALLBACK(void) callback(const char*));
    C_API_EXPORT(void) carl_tickCallbacks(uint64_t sessionPtr);
    C_API_EXPORT(void) carl_addSerializedInputSample(uint64_t sessionPtr, uint8_t* bytes, uint64_t size);
    C_API_EXPORT(void) carl_addInputSample(uint64_t sessionPtr, carl_InputSample* sample);
    C_API_EXPORT(void) carl_disposeSession(uint64_t sessionPtr);
    C_API_EXPORT(void) carl_createRecognizerAsync(uint64_t sessionPtr, uint64_t definitionPtr, uint64_t requestId, C_API_CALLBACK(void) callback(uint64_t, uint64_t));
    C_API_EXPORT(double) carl_getCurrentScore(uint64_t recognizerPtr);
    C_API_EXPORT(void) carl_setSensitivity(uint64_t recognizerPtr, double sensitivity);
    C_API_EXPORT(uint64_t) carl_getCanonicalRecordingInspector(uint64_t recognizerPtr);
    C_API_EXPORT(void) carl_disposeRecognizer(uint64_t sessionPtr, uint64_t recognizerPtr);
}