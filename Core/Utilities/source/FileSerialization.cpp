/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "carl/utilities/FileSerialization.h"

#include <carl/Definition.h>
#include <carl/Example.h>
#include <carl/Serialization.h>

#include <fstream>

namespace carl::utilities
{
    namespace
    {
        constexpr uint32_t VERSION{ 1 };

        struct FileHeader
        {
            enum class FileType : uint32_t
            {
                Undefined = 0,
                Example = 1,
                Definition = 2,
                COUNT,
            };

            std::array<char, 4> Magic{ 'C', 'A', 'R', 'L' };
            uint32_t Version{ VERSION };
            FileType Type{ FileType::Undefined };

            bool isValid(FileType type) const
            {
                return
                    Magic[0] == 'C' &&
                    Magic[1] == 'A' &&
                    Magic[2] == 'R' &&
                    Magic[3] == 'L' &&
                    Version == VERSION &&
                    Type == type;
            }
        };

        /// <summary>
        /// The very first implementation of InputSample, from commit 6d77871.
        /// These input legacySamples are left-handed, using Unity conventions, and
        /// rely on a separate joint set based on the old OVR joints.
        /// </summary>
        struct LegacyInputSampleV0
        {
            enum class Joint : uint64_t
            {
                UNUSED_HandJointId_HandThumb0,
                ThumbFingerBase,
                UNUSED_HandJointId_HandThumb2,
                UNUSED_HandJointId_HandThumb3,
                ThumbFingerTip,
                IndexFingerBase,
                UNUSED_HandJointId_HandIndex2,
                UNUSED_HandJointId_HandIndex3,
                IndexFingerTip,
                MiddleFingerBase,
                UNUSED_HandJointId_HandMiddle2,
                UNUSED_HandJointId_HandMiddle3,
                MiddleFingerTip,
                RingFingerBase,
                UNUSED_HandJointId_HandRing2,
                UNUSED_HandJointId_HandRing3,
                RingFingerTip,
                UNUSED_HandJointId_HandPinky0,
                LittleFingerBase,
                UNUSED_HandJointId_HandPinky2,
                UNUSED_HandJointId_HandPinky3,
                LittleFingerTip,
                COUNT
            };

            static constexpr std::array<int, 26> LEGACY_TO_MODERN_HAND_JOINT_MAPPING
            {
                0,
                -1,
                1,
                2,
                3,
                4,
                -1,
                5,
                6,
                7,
                8,
                -1,
                9,
                10,
                11,
                12,
                -1,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
                21,
            };

            LegacyInputSampleV0() = default;
            LegacyInputSampleV0(Deserialization& deserialization)
            {
                deserialization >> Timestamp;
                deserialization >> HmdPose;
                deserialization >> LeftWristPose;
                deserialization >> RightWristPose;
                deserialization >> LeftHandJointPoses;
                deserialization >> RightHandJointPoses;
            }

            operator InputSample() const
            {
                constexpr auto switchHandedness = [](const TransformT& leftHanded) {
                    QuaternionT rotation{ leftHanded.rotation() };
                    rotation.z() *= -1;
                    rotation.w() *= -1;
                    VectorT translation{ leftHanded.translation() };
                    translation.z() *= -1;
                    TransformT rightHanded{};
                    rightHanded.fromPositionOrientationScale(translation, rotation, UNIT_SCALE);
                    return rightHanded;
                };

                InputSample sample{};
                sample.Timestamp = Timestamp;
                if (HmdPose.has_value())
                {
                    sample.HmdPose = switchHandedness(HmdPose.value());
                }
                if (LeftWristPose.has_value())
                {
                    sample.LeftWristPose = switchHandedness(LeftWristPose.value());
                }
                if (RightWristPose.has_value())
                {
                    sample.RightWristPose = switchHandedness(RightWristPose.value());
                }
                if (LeftHandJointPoses.has_value())
                {
                    sample.LeftHandJointPoses.emplace();
                    for (size_t idx = 0; idx < LEGACY_TO_MODERN_HAND_JOINT_MAPPING.size(); ++idx)
                    {
                        if (LEGACY_TO_MODERN_HAND_JOINT_MAPPING[idx] >= 0)
                        {
                            sample.LeftHandJointPoses.value()[idx] = switchHandedness(LeftHandJointPoses.value()[LEGACY_TO_MODERN_HAND_JOINT_MAPPING[idx]]);
                        }
                        else
                        {
                            sample.LeftHandJointPoses.value()[idx] = TransformT::Identity();
                        }
                    }
                }
                if (HmdPose.has_value())
                {
                    sample.RightHandJointPoses.emplace();
                    for (size_t idx = 0; idx < LEGACY_TO_MODERN_HAND_JOINT_MAPPING.size(); ++idx)
                    {
                        if (LEGACY_TO_MODERN_HAND_JOINT_MAPPING[idx] >= 0)
                        {
                            sample.RightHandJointPoses.value()[idx] = switchHandedness(RightHandJointPoses.value()[LEGACY_TO_MODERN_HAND_JOINT_MAPPING[idx]]);
                        }
                        else
                        {
                            sample.RightHandJointPoses.value()[idx] = TransformT::Identity();
                        }
                    }
                }
                return sample;
            }

            double Timestamp{};
            std::optional<TransformT> HmdPose{};
            std::optional<TransformT> LeftWristPose{};
            std::optional<TransformT> RightWristPose{};
            std::optional<std::array<TransformT, static_cast<size_t>(Joint::COUNT)>> LeftHandJointPoses{};
            std::optional<std::array<TransformT, static_cast<size_t>(Joint::COUNT)>> RightHandJointPoses{};
        };

        /// <summary>
        /// When InputSample switched to OpenXR conventions in 1c47934.
        /// These input legacySamples are right-handed, using OpenXR conventions.
        /// </summary>
        struct LegacyInputSampleV1
        {
            LegacyInputSampleV1() = default;
            LegacyInputSampleV1(Deserialization& deserialization)
            {
                deserialization >> Timestamp;
                deserialization >> HmdPose;
                deserialization >> LeftWristPose;
                deserialization >> RightWristPose;
                deserialization >> LeftHandJointPoses;
                deserialization >> RightHandJointPoses;
            }

            operator InputSample() const
            {
                InputSample sample{};
                sample.Timestamp = Timestamp;
                sample.HmdPose = HmdPose;
                sample.LeftWristPose = LeftWristPose;
                sample.RightWristPose = RightWristPose;
                sample.LeftHandJointPoses = LeftHandJointPoses;
                sample.RightHandJointPoses = RightHandJointPoses;

                return sample;
            }

            double Timestamp{};
            std::optional<TransformT> HmdPose{};
            std::optional<TransformT> LeftWristPose{};
            std::optional<TransformT> RightWristPose{};
            std::optional<std::array<TransformT, static_cast<size_t>(InputSample::Joint::COUNT)>> LeftHandJointPoses{};
            std::optional<std::array<TransformT, static_cast<size_t>(InputSample::Joint::COUNT)>> RightHandJointPoses{};
        };

        std::optional<action::InProgressRecording> TryDeserializeLegacyRecording(Deserialization& deserialization)
        {
            uint64_t recordingSize{};
            deserialization >> recordingSize;

            std::optional<action::InProgressRecording> recording{};

            if (!recording.has_value())
            {
                try
                {
                    carl::Deserialization d{ deserialization };
                    std::vector<LegacyInputSampleV0> legacySamples{};
                    legacySamples.reserve(recordingSize);
                    for (size_t idx = 0; idx < recordingSize; ++idx)
                    {
                        legacySamples.emplace_back(d);
                    }
                    bool samplesValid = true;
                    for (int idx = 1; samplesValid && idx < legacySamples.size(); ++idx)
                    {
                        auto delta = legacySamples[idx].Timestamp - legacySamples[idx - 1].Timestamp;
                        samplesValid = delta > 0 && delta < 60;
                    }
                    if (samplesValid)
                    {
                        recording.emplace();
                        for (const auto& legacySample : legacySamples)
                        {
                            recording->addSample(legacySample);
                        }
                        deserialization = d;
                    }
                }
                catch (...) {}
            }

            if (!recording.has_value())
            {
                try
                {
                    carl::Deserialization d{ deserialization };
                    std::vector<LegacyInputSampleV1> legacySamples{};
                    legacySamples.reserve(recordingSize);
                    for (size_t idx = 0; idx < recordingSize; ++idx)
                    {
                        legacySamples.emplace_back(d);
                    }
                    bool samplesValid = true;
                    for (int idx = 1; samplesValid && idx < legacySamples.size(); ++idx)
                    {
                        auto delta = legacySamples[idx].Timestamp - legacySamples[idx - 1].Timestamp;
                        samplesValid = delta > 0 && delta < 60;
                    }
                    if (samplesValid)
                    {
                        recording.emplace();
                        for (const auto& legacySample : legacySamples)
                        {
                            recording->addSample(legacySample);
                        }
                        deserialization = d;
                    }
                }
                catch (...) {}
            }

            if (!recording.has_value())
            {
                try
                {
                    carl::Deserialization d{ deserialization };
                    std::vector<InputSample> legacySamples{};
                    legacySamples.reserve(recordingSize);
                    for (size_t idx = 0; idx < recordingSize; ++idx)
                    {
                        legacySamples.emplace_back(d);
                    }
                    bool samplesValid = true;
                    for (int idx = 1; samplesValid && idx < legacySamples.size(); ++idx)
                    {
                        auto delta = legacySamples[idx].Timestamp - legacySamples[idx - 1].Timestamp;
                        samplesValid = delta > 0 && delta < 60;
                    }
                    if (samplesValid)
                    {
                        recording.emplace();
                        for (const auto& legacySample : legacySamples)
                        {
                            recording->addSample(legacySample);
                        }
                        deserialization = d;
                    }
                }
                catch (...) {}
            }

            return recording;
        }

        std::vector<uint8_t> ReadFileBytes(std::filesystem::path path)
        {
            auto length = std::filesystem::file_size(path);

            std::ifstream fileStream{ path.c_str(), std::ios::in | std::ios::binary };
            std::vector<uint8_t> bytes{};
            bytes.resize(static_cast<size_t>(length));
            fileStream.seekg(std::ios::beg);
            fileStream.read(reinterpret_cast<char*>(bytes.data()), bytes.size());

            return bytes;
        }
    }

    template<>
    std::vector<uint8_t> Serialize(const carl::action::Example& example)
    {
        FileHeader header{};
        header.Type = FileHeader::FileType::Example;

        std::vector<uint8_t> bytes{};
        carl::Serialization serialization{ bytes };

        serialization << header;
        example.serialize(serialization);

        return bytes;
    }

    template<>
    void SerializeToFile(const carl::action::Example& example, std::filesystem::path path)
    {
        std::vector<uint8_t> bytes = Serialize(example);

        std::ofstream fileStream{ path, std::ios::out | std::ios::binary };
        fileStream.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
        fileStream.close();
    }

    template<>
    std::optional<carl::action::Example> TryDeserialize(gsl::span<const uint8_t> bytes)
    {
        carl::Deserialization deserialization{ bytes };

        FileHeader header{};
        deserialization >> header;
        if (!header.isValid(FileHeader::FileType::Example))
        {
            return TryDeserializeLegacy<action::Example>(bytes);
        }

        return { { deserialization } };
    }

    template<>
    std::optional<carl::action::Example> TryDeserializeFromFile(std::filesystem::path path)
    {
        std::vector<uint8_t> bytes = ReadFileBytes(path);
        return TryDeserialize<action::Example>(bytes);
    }

    template<>
    std::optional<carl::action::Example> TryDeserializeLegacy(gsl::span<const uint8_t> bytes)
    {
        carl::Deserialization deserialization{ bytes };
        float startTimestamp{};
        deserialization >> startTimestamp;
        float endTimestamp{};
        deserialization >> endTimestamp;
        int recordingBytesLength{};
        deserialization >> recordingBytesLength;

        float duration = endTimestamp - startTimestamp;
        if (duration < 0.f || duration > 1000.f)
        {
            // Early-out for insane timestamps.
            return {};
        }

        if (recordingBytesLength != static_cast<int>(bytes.size()) - 12)
        {
            // Early-out if recording bytes length is wrong, which indicates that this is not a supported legacy example file.
            return {};
        }

        auto recording = TryDeserializeLegacyRecording(deserialization);

        if (!recording.has_value() || deserialization.remainingBytes() != 0)
        {
            return{};
        }
        return{ action::Example{ { std::move(recording.value()) }, startTimestamp, endTimestamp } };
    }

    template<>
    std::optional<carl::action::Example> TryDeserializeLegacyFile(std::filesystem::path path)
    {
        std::vector<uint8_t> bytes = ReadFileBytes(path);
        return TryDeserializeLegacy<action::Example>(bytes);
    }

    template<>
    std::vector<uint8_t> Serialize(const carl::action::Definition& definition)
    {
        FileHeader header{};
        header.Type = FileHeader::FileType::Definition;

        std::vector<uint8_t> bytes{};
        carl::Serialization serialization{ bytes };

        serialization << header;
        definition.serialize(serialization);

        return bytes;
    }

    template<>
    void SerializeToFile(const carl::action::Definition& definition, std::filesystem::path path)
    {
        std::vector<uint8_t> bytes = Serialize(definition);

        std::ofstream fileStream{ path, std::ios::out | std::ios::binary };
        fileStream.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
        fileStream.close();
    }

    template<>
    std::optional<carl::action::Definition> TryDeserialize(gsl::span<const uint8_t> bytes)
    {
        carl::Deserialization deserialization{ bytes };

        FileHeader header{};
        deserialization >> header;
        if (!header.isValid(FileHeader::FileType::Definition))
        {
            return TryDeserializeLegacy<action::Definition>(bytes);
        }

        return { { deserialization } };
    }

    template<>
    std::optional<carl::action::Definition> TryDeserializeFromFile(std::filesystem::path path)
    {
        std::vector<uint8_t> bytes = ReadFileBytes(path);
        return TryDeserialize<action::Definition>(bytes);
    }

    template<>
    std::optional<carl::action::Definition> TryDeserializeLegacy(gsl::span<const uint8_t> bytes)
    {
        carl::Deserialization deserialization{ bytes };

        int fileSize{};
        deserialization >> fileSize;

        if (fileSize != static_cast<int>(bytes.size()) - 4)
        {
            return {};
        }

        action::Definition::ActionType actionType{};
        deserialization >> actionType;
        action::Definition definition{ actionType };

        uint64_t count{};
        deserialization >> count;
        for (uint64_t idx = 0; idx < count; ++idx)
        {
            auto recording = TryDeserializeLegacyRecording(deserialization);
            if (!recording.has_value())
            {
                return{};
            }
            double startTimestamp{};
            deserialization >> startTimestamp;
            double endTimestamp{};
            deserialization >> endTimestamp;
            definition.addExample({ std::move(recording.value()), startTimestamp, endTimestamp });
        }
        deserialization >> count;
        std::vector<action::Example> counterexamples{};
        counterexamples.reserve(count);
        for (uint64_t idx = 0; idx < count; ++idx)
        {
            auto recording = TryDeserializeLegacyRecording(deserialization);
            if (!recording.has_value())
            {
                return{};
            }
            double startTimestamp{};
            deserialization >> startTimestamp;
            double endTimestamp{};
            deserialization >> endTimestamp;
            definition.addCounterexample({ std::move(recording.value()), startTimestamp, endTimestamp });
        }
        deserialization >> definition.DefaultSensitivity;

        if (deserialization.remainingBytes() != 0)
        {
            return{};
        }
        return definition;
    }

    template<>
    std::optional<carl::action::Definition> TryDeserializeLegacyFile(std::filesystem::path path)
    {
        std::vector<uint8_t> bytes = ReadFileBytes(path);
        return TryDeserializeLegacy<action::Definition>(bytes);
    }
}
