/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/InputSample.h>

#include <gsl/span>

#include <cassert>

namespace carl
{
    namespace
    {
        template<typename T, size_t Count>
        void lerpOptionalArray(
            const std::optional<std::array<T, Count>>& l, 
            const std::optional<std::array<T, Count>>& r, 
            double t, 
            std::optional<std::array<T, Count>>& target)
        {
            if (l.has_value() && r.has_value())
            {
                gsl::span<const T> lSamples{ l.value() };
                gsl::span<const T> rSamples{ r.value() };
                assert(lSamples.size() == rSamples.size());
                std::array<T, Count> lerped{};
                for (size_t idx = 0; idx < lSamples.size(); ++idx)
                {
                    lerped[idx] = math::Lerp(lSamples[idx], rSamples[idx], static_cast<NumberT>(t));
                }
                target = std::move(lerped);
            }
            else if (l.has_value())
            {
                target = l.value();
            }
            else if (r.has_value())
            {
                target = r.value();
            }
        }
    }

    InputSample::InputSample(Deserialization& deserialization)
    {
        deserialization >> Timestamp;
        deserialization >> HmdPose;
        deserialization >> LeftWristPose;
        deserialization >> RightWristPose;
        deserialization >> LeftHandJointPoses;
        deserialization >> RightHandJointPoses;
        deserialization >> LeftControllerInput;
        deserialization >> RightControllerInput;
    }

    void InputSample::serialize(Serialization& serialization) const
    {
        serialization << Timestamp;
        serialization << HmdPose;
        serialization << LeftWristPose;
        serialization << RightWristPose;
        serialization << LeftHandJointPoses;
        serialization << RightHandJointPoses;
        serialization << LeftControllerInput;
        serialization << RightControllerInput;
    }

    InputSample InputSample::Lerp(const InputSample& a, const InputSample& b, double t)
    {
        InputSample result{};

        result.Timestamp = (1. - t) * a.Timestamp + t * b.Timestamp;

        constexpr auto lerpOptional =
            [](const auto& l, const auto& r, double t) -> std::optional<TransformT> {
            if (l.has_value() && r.has_value())
            {
                return math::Lerp(l.value(), r.value(), static_cast<NumberT>(t));
            }
            else if (l.has_value())
            {
                return l;
            }
            else
            {
                return r;
            }
        };

        result.HmdPose = lerpOptional(a.HmdPose, b.HmdPose, t);
        result.LeftWristPose = lerpOptional(a.LeftWristPose, b.LeftWristPose, t);
        result.RightWristPose = lerpOptional(a.RightWristPose, b.RightWristPose, t);

        lerpOptionalArray(a.LeftHandJointPoses, b.LeftHandJointPoses, t, result.LeftHandJointPoses);
        lerpOptionalArray(a.RightHandJointPoses, b.RightHandJointPoses, t, result.RightHandJointPoses);

        lerpOptionalArray(a.LeftControllerInput, b.LeftControllerInput, t, result.LeftControllerInput);
        lerpOptionalArray(a.RightControllerInput, b.RightControllerInput, t, result.RightControllerInput);

        return result;
    }
}
