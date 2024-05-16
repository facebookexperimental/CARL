/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "EgocentricTemporalSpace.h"

#include "carl/Example.h"
#include "DynamicTimeWarping.h"

#include <gsl/span>

#include <array>

namespace carl::descriptor
{
    // This namespace contains small, self-contained math constructs for use ONLY within
    // descriptors. These are needed because Eigen matrices are not trivially copyable,
    // so in order for descriptors to remain trivially copyable they need math primitives
    // that are also trivially copyable. Note that these types should be kept as minimal
    // as possible, containing only what's needed for descriptor comparison; Eigen should
    // still be used for all other math.
    namespace trivial
    {
        struct Point
        {
            NumberT X{};
            NumberT Y{};
            NumberT Z{};

            Point() = default;

            Point(const VectorT& v)
                : X{ v.x() }
                , Y{ v.y() }
                , Z{ v.z() }
            {
            }

            Point(NumberT x, NumberT y, NumberT z)
                : X{ x }
                , Y{ y }
                , Z{ z }
            {
            }

            Point& operator=(const VectorT& v)
            {
                X = v.x();
                Y = v.y();
                Z = v.z();
                return *this;
            }

            operator VectorT() const
            {
                return{ X, Y, Z };
            }

            Point operator-(const Point& other) const
            {
                return{ X - other.X, Y - other.Y, Z - other.Z };
            }

            NumberT lengthSq() const
            {
                return X * X + Y * Y + Z * Z;
            }

            NumberT length() const
            {
                return std::sqrt(lengthSq());
            }

            NumberT distanceSq(const Point& other) const
            {
                return (operator-(other)).lengthSq();
            }

            NumberT distance(const Point& other) const
            {
                return std::sqrt(distanceSq(other));
            }
        };

        struct Quaternion
        {
            NumberT X{};
            NumberT Y{};
            NumberT Z{};
            NumberT W{};

            Quaternion() = default;

            Quaternion(const QuaternionT& q)
                : X{ q.x() }
                , Y{ q.y() }
                , Z{ q.z() }
                , W{ q.w() }
            {
            }

            Quaternion(NumberT x, NumberT y, NumberT z, NumberT w)
                : X{ x }
                , Y{ y }
                , Z{ z }
                , W{ w }
            {
            }

            Quaternion& operator=(const QuaternionT& other)
            {
                X = other.x();
                Y = other.y();
                Z = other.z();
                W = other.w();
                return *this;
            }

            Quaternion operator*(const Quaternion& other) const
            {
                return{
                    W * other.X + X * other.W + Y * other.Z - Z * other.Y,
                    W * other.Y - X * other.Z + Y * other.W + Z * other.X,
                    W * other.Z + X * other.Y - Y * other.X + Z * other.W,
                    W * other.W - X * other.X - Y * other.Y - Z * other.Z,
                };
            }

            operator QuaternionT() const
            {
                return{ W, X, Y, Z };
            }

            Quaternion conjugate() const
            {
                return{ -X, -Y, -Z, W };
            }

            NumberT angularDistance(const Quaternion& other) const
            {
                auto qt{ operator*(other.conjugate()) };
                Point pt{ qt.X, qt.Y, qt.Z };
                return static_cast<NumberT>(2) * std::atan2(pt.length(), std::abs(qt.W));
            }
        };
    }

    template<typename...>
    struct ArraySizeCount;

    template<typename T, size_t S>
    struct ArraySizeCount<std::array<T, S>>
    {
        static constexpr auto Size = S;
    };

    template<typename T, size_t S, typename... Ts>
    struct ArraySizeCount<std::array<T, S>, Ts...>
    {
        static constexpr auto Size = S + ArraySizeCount<Ts...>::Size;
    };

    template<typename T, size_t l, size_t r>
    constexpr std::array<T, l + r> arrayConcat(const std::array<T, l>& a, const std::array<T, r>& b)
    {
        std::array<T, l + r> ret{};
        for (int idx = 0; idx < l; ++idx)
        {
            ret[idx] = a[idx];
        }

        for (int idx = 0; idx < r; ++idx)
        {
            ret[idx + l] = b[idx];
        }

        return ret;
    }

    template<typename T, size_t l, size_t r, typename... Ts>
    constexpr std::array<T, l + r + ArraySizeCount<Ts...>::Size> arrayConcat(
        const std::array<T, l>& a,
        const std::array<T, r>& b,
        const Ts&... ts)
    {
        return arrayConcat(a, arrayConcat(b, ts...));
    }

    template<typename...>
    struct Tuning;

    template<typename T1, typename T2>
    struct Tuning<T1, T2>
    {
        static constexpr auto DEFAULT_TUNING{ arrayConcat(T1::DEFAULT_TUNING, T2::DEFAULT_TUNING) };

        template<typename T>
        static gsl::span<const NumberT> getTuning(gsl::span<const NumberT> values)
        {
            if constexpr (std::is_same<T, T1>::value)
            {
                return gsl::make_span(values.data(), T1::DEFAULT_TUNING.size());
            }
            else if constexpr (std::is_same<T, T2>::value)
            {
                return gsl::make_span(values.data() + T1::DEFAULT_TUNING.size(), T2::DEFAULT_TUNING.size());
            }
        }
    };

    template<typename T1, typename T2, typename... Ts>
    struct Tuning<T1, T2, Ts...>
    {
        static constexpr auto DEFAULT_TUNING{ arrayConcat(T1::DEFAULT_TUNING, Tuning<T2, Ts...>::DEFAULT_TUNING) };

        template<typename T>
        static gsl::span<const NumberT> getTuning(gsl::span<const NumberT> values)
        {
            constexpr auto size = T1::DEFAULT_TUNING.size();
            if constexpr (std::is_same<T, T1>::value)
            {
                return gsl::make_span(values.data(), size);
            }
            else
            {
                return Tuning<T2, Ts...>::template getTuning<T>(gsl::make_span(values.data() + size, values.size() - size));
            }
        }
    };

    enum class Handedness : uint64_t
    {
        LeftHanded,
        RightHanded,
        TwoHanded,
    };

    template<Handedness Handedness>
    inline TransformT getWristPose(const InputSample& sample)
    {
        constexpr size_t Y_AXIS_JOINT{
            static_cast<size_t>(InputSample::Joint::LittleFingerBase) };
        constexpr size_t Z_AXIS_JOINT{
            static_cast<size_t>(InputSample::Joint::IndexFingerBase) };

        constexpr bool isLeftHanded = Handedness == Handedness::LeftHanded;
        const auto& wristPose =
            isLeftHanded ? sample.LeftWristPose.value() : sample.RightWristPose.value();
        const auto& jointPoses =
            isLeftHanded ? sample.LeftHandJointPoses.value() : sample.RightHandJointPoses.value();

        auto zAxis = (jointPoses[Z_AXIS_JOINT].translation() - wristPose.translation()).normalized();
        auto yAxis = zAxis.cross((jointPoses[Y_AXIS_JOINT].translation() - wristPose.translation()).cross(zAxis)).normalized();
        return math::LookTransform(zAxis, yAxis, wristPose.translation());
    }

    //#define DESCRIPTOR_ANALYSIS
#ifdef DESCRIPTOR_ANALYSIS
    constexpr NumberT NULL_TUNING{ 0 };
    constexpr auto createDistanceFunction(NumberT identicalityThreshold, NumberT irreconcilabilityThreshold)
    {
        return [identicalityThreshold, irreconcilabilityThreshold](NumberT distance, NumberT tuning)
            {
                return tuning * distance;
            };
    }
#else
    constexpr NumberT NULL_TUNING{ 1000000 };
    constexpr auto createDistanceNormalizationFunction(NumberT identicalityThreshold, NumberT irreconcilabilityThreshold)
    {
        return [identicalityThreshold, irreconcilabilityThreshold](NumberT distance, NumberT tuning)
            {
                auto lowerBound = tuning * identicalityThreshold;
                auto upperBound = tuning * irreconcilabilityThreshold;
                return std::max<NumberT>(0, std::pow<NumberT>((distance - lowerBound) / (upperBound - lowerBound), 3));
            };
    }
#endif

    // TODO: Find a better place for this (and everything above it) to live.
    template<typename DescriptorT>
    void extendSequence(const InputSample& newSample, std::vector<DescriptorT>& sequence, InputSample& mostRecentSample, gsl::span<const NumberT> tuning)
    {
        constexpr NumberT THRESHOLD{ 0.5 };

        // Handle startup, in which case m_mostRecentSample will be a default value.
        if (sequence.empty())
        {
            auto descriptor = DescriptorT::TryCreate(newSample, newSample);
            if (descriptor.has_value())
            {
                sequence.emplace_back(std::move(*descriptor));
                mostRecentSample = newSample;
            }
        }
        else
        {
            // Iteratively create new mostRecentSamples and descriptors until the most recent descriptor is sufficiently close to that of newSample.
            // TODO: Figure out if delta descriptors (rotation, translation, etc.) will play correctly with this, since lerping inputs should not 
            //       change them. Might be necessary to mute such descriptors with tuning.
            while (true)
            {
                auto sampleDesc = DescriptorT::TryCreate(newSample, mostRecentSample);

                // If we weren't able to create a descriptor, stop iterating.
                if (!sampleDesc.has_value())
                {
                    break;
                }

                auto outerDistance = DescriptorT::Distance(*sampleDesc, sequence.back(), tuning);
                if (outerDistance < THRESHOLD)
                {
                    break;
                }

                // Posit new samples, and associated descriptors, to bring the end of sequence closer to newSample.
                NumberT upper = 1;
                NumberT lower = 0;
                NumberT mid{};
                while (true)
                {
                    mid = (upper + lower) / NumberT{ 2 };
                    auto intermediateDesc = DescriptorT::Lerp(sequence.back(), *sampleDesc, mid);
                    auto distance = DescriptorT::Distance(intermediateDesc, sequence.back(), tuning);
                    if (distance > THRESHOLD)
                    {
                        // intermediate sample is too distant, continue searching for a nearer sample
                        upper = mid;
                    }
                    else if (distance < NumberT{ 0.9 } *THRESHOLD)
                    {
                        // intermediate sample is too close, continue searching for a more distant sample
                        lower = mid;
                    }
                    else
                    {
                        // intermediate sample is satisfactory, stop searching
                        sequence.push_back(intermediateDesc);
                        mostRecentSample = InputSample::Lerp(mostRecentSample, newSample, mid);
                        break;
                    }
                }
            }
        }
    }

    template<typename DescriptorT>
    std::vector<DescriptorT> createDescriptorSequenceFromRecording(const action::Recording& recording, double startTimestamp, double endTimestamp, gsl::span<const NumberT> tuning)
    {
        std::vector<DescriptorT> sequence{};

        auto samples = recording.getSamples();
        InputSample mostRecentSample{};

        size_t idx = 0;
        while (idx < samples.size() - 1 && samples[idx + 1].Timestamp < startTimestamp)
        {
            ++idx;
        }
        while (idx < samples.size() - 1 && samples[idx + 1].Timestamp < endTimestamp)
        {
            descriptor::extendSequence(samples[idx], sequence, mostRecentSample, tuning);
            ++idx;
        }

        return sequence;
    }

    template<typename DescriptorT>
    std::vector<std::vector<DescriptorT>> createDescriptorSequencesFromExamples(gsl::span<const action::Example> examples, gsl::span<const NumberT> tuning, double padding = 0.)
    {
        std::vector<std::vector<DescriptorT>> sequences{};
        sequences.reserve(examples.size());
        for (const auto& example : examples)
        {
            sequences.push_back(createDescriptorSequenceFromRecording<DescriptorT>(example.getRecording(), example.getStartTimestamp() - padding, example.getEndTimestamp() + padding, tuning));
        }
        return sequences;
    }

    template<Handedness Handedness>
    class HandShape
    {
        static constexpr std::array<size_t, 5> JOINTS{
            static_cast<size_t>(InputSample::Joint::ThumbFingerTip),
                static_cast<size_t>(InputSample::Joint::IndexFingerTip),
                static_cast<size_t>(InputSample::Joint::MiddleFingerTip),
                static_cast<size_t>(InputSample::Joint::RingFingerTip),
                static_cast<size_t>(InputSample::Joint::LittleFingerTip),
        };
        static constexpr size_t NORMALIZATION_JOINT{
            static_cast<size_t>(InputSample::Joint::IndexFingerBase) };

    public:
        static constexpr std::array<NumberT, JOINTS.size()> DEFAULT_TUNING{ 1., 1., 1., 1., 1. };
        static constexpr auto HANDEDNESS{ Handedness };

        static std::optional<HandShape> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if constexpr (Handedness == Handedness::LeftHanded)
            {
                if (sample.LeftWristPose.has_value() && sample.LeftHandJointPoses.has_value())
                {
                    return HandShape{ sample, priorSample };
                }
            }
            else if constexpr (Handedness == Handedness::RightHanded) {
                if (sample.RightWristPose.has_value() && sample.RightHandJointPoses.has_value()) {
                    return HandShape{ sample, priorSample };
                }
            }
            return {};
        }

        static NumberT Distance(const HandShape& a, const HandShape& b, gsl::span<const NumberT> tuning)
        {
            NumberT distance = 0;
            for (size_t idx = 0; idx < a.m_positions.size(); ++idx) {
                distance = std::max(distance, normalizeDistance(a.m_positions[idx].distance(b.m_positions[idx]), tuning[idx]));
            }
            return distance;
        }

        static HandShape Lerp(const HandShape& a, const HandShape& b, NumberT t)
        {
            HandShape result{};
            for (size_t idx = 0; idx < JOINTS.size(); ++idx)
            {
                VectorT aVec{ a.m_positions[idx] };
                VectorT bVec{ b.m_positions[idx] };
                result.m_positions[idx] = (static_cast<NumberT>(1) - t) * aVec + t * bVec;
            }
            return result;
        }

        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(gsl::span<const action::Example> examples)
        {
            if (examples.size() < 2)
            {
                return DEFAULT_TUNING;
            }

            auto sequences = createDescriptorSequencesFromExamples<HandShape<Handedness>>(examples, DEFAULT_TUNING);
            auto extendedSequences = createDescriptorSequencesFromExamples<HandShape<Handedness>>(examples, DEFAULT_TUNING, 1.);

            auto tuning = DEFAULT_TUNING;
            for (size_t idx = 0; idx < tuning.size(); ++idx)
            {
                NumberT maxConnectionCost = 0;
                for (const auto& extendedSequence : extendedSequences)
                {
                    for (const auto& sequence : sequences)
                    {
                        auto distanceFunction = [idx](const auto& a, const auto& b) {
                            return a.m_positions[idx].distance(b.m_positions[idx]);
                        };
                        auto result = DynamicTimeWarping::Match<const HandShape<Handedness>>(extendedSequence, sequence, distanceFunction);
                        maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                    }
                }
                NumberT midpoint = (IDENTICALITY_THRESHOLD + IRRECONCILABILITY_THRESHOLD) / 2;
                tuning[idx] = std::max(maxConnectionCost / midpoint, DEFAULT_TUNING[idx]);
            }
            return tuning;
        }

        HandShape() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.01 };
        static inline constexpr NumberT IRRECONCILABILITY_THRESHOLD{ 0.03 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD, IRRECONCILABILITY_THRESHOLD) };
        static inline constexpr NumberT CANONICAL_NORMALIZATION_LENGTH{ 0.1 };
        std::array<trivial::Point, JOINTS.size()> m_positions{};

        HandShape(const InputSample& sample, const InputSample&)
        {
            constexpr bool isLeftHanded = Handedness == Handedness::LeftHanded;
            const auto& jointPoses =
                isLeftHanded ? sample.LeftHandJointPoses.value() : sample.RightHandJointPoses.value();

            auto inverseWristPose = getWristPose<Handedness>(sample).inverse();

            NumberT normalizationFactor = CANONICAL_NORMALIZATION_LENGTH /
                (inverseWristPose * jointPoses[NORMALIZATION_JOINT].translation()).norm();
            for (size_t idx = 0; idx < JOINTS.size(); ++idx)
            {
                m_positions[idx] = normalizationFactor * (inverseWristPose * jointPoses[JOINTS[idx]].translation());
            }
        }
    };

    template<Handedness Handedness>
    class EgocentricWristOrientation
    {
    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };
        static constexpr auto HANDEDNESS{ Handedness };

        static std::optional<EgocentricWristOrientation> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if (!sample.HmdPose.has_value())
            {
                return{};
            }

            if constexpr (Handedness == Handedness::LeftHanded)
            {
                if (sample.LeftWristPose.has_value() && priorSample.LeftWristPose.has_value())
                {
                    return EgocentricWristOrientation{ sample, priorSample };
                }
            }
            else if constexpr (Handedness == Handedness::RightHanded)
            {
                if (sample.RightWristPose.has_value() && priorSample.RightWristPose.has_value())
                {
                    return EgocentricWristOrientation{ sample, priorSample };
                }
            }
            return{};
        }

        static NumberT Distance(
            const EgocentricWristOrientation& a,
            const EgocentricWristOrientation& b,
            gsl::span<const NumberT> tuning)
        {
            auto distance = a.m_egocentricTemporalOrientation.angularDistance(b.m_egocentricTemporalOrientation);
            return normalizeDistance(distance, tuning[0]);
        }

        static EgocentricWristOrientation Lerp(
            const EgocentricWristOrientation& a,
            const EgocentricWristOrientation& b,
            NumberT t)
        {
            QuaternionT aQuat{ a.m_egocentricTemporalOrientation };
            QuaternionT bQuat{ b.m_egocentricTemporalOrientation };
            EgocentricWristOrientation result{};
            result.m_egocentricTemporalOrientation = aQuat.slerp(t, bQuat);
            return result;
        }

        template<typename ExamplesT>
        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(const ExamplesT& examples)
        {
            if (examples.size() < 2)
            {
                return DEFAULT_TUNING;
            }

            auto sequences = createDescriptorSequencesFromExamples<EgocentricWristOrientation<Handedness>>(examples, DEFAULT_TUNING);
            auto extendedSequences = createDescriptorSequencesFromExamples<EgocentricWristOrientation<Handedness>>(examples, DEFAULT_TUNING, 1.);

            auto tuning = DEFAULT_TUNING;
            NumberT maxConnectionCost = 0;
            for (const auto& extendedSequence : extendedSequences)
            {
                for (const auto& sequence : sequences)
                {
                    auto distanceFunction = [](const auto& a, const auto& b) {
                        return a.m_egocentricTemporalOrientation.angularDistance(b.m_egocentricTemporalOrientation);
                    };
                    auto result = DynamicTimeWarping::Match<const EgocentricWristOrientation<Handedness>>(extendedSequence, sequence, distanceFunction);
                    maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                }
            }
            NumberT midpoint = (IDENTICALITY_THRESHOLD + IRRECONCILABILITY_THRESHOLD) / 2;
            tuning[0] = std::max(maxConnectionCost / midpoint, DEFAULT_TUNING[0]);
            return tuning;
        }

        EgocentricWristOrientation() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.17453 };
        static inline constexpr NumberT IRRECONCILABILITY_THRESHOLD{ 0.5236 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD, IRRECONCILABILITY_THRESHOLD) };
        trivial::Quaternion m_egocentricTemporalOrientation{};

        EgocentricWristOrientation(const InputSample& sample, const InputSample&)
        {
            constexpr bool isLeftHanded = Handedness == Handedness::LeftHanded;
            auto wristPose = getWristPose<Handedness>(sample);
            auto wristPosition = wristPose.translation();
            auto& hmdPose = sample.HmdPose.value();

            auto ets = EgocentricTemporalSpace::getPose(wristPosition, hmdPose);
            m_egocentricTemporalOrientation = QuaternionT{ ets.rotation().inverse() * wristPose.rotation() };
        }
    };

    template<Handedness Handedness>
    class HandPose
    {
    public:
        using TuningT = Tuning<HandShape<Handedness>, EgocentricWristOrientation<Handedness>>;
        static constexpr auto DEFAULT_TUNING{ TuningT::DEFAULT_TUNING };
        static constexpr auto HANDEDNESS{ Handedness };

        static std::optional<HandPose> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            auto handShapeSample = HandShape<Handedness>::TryCreate(sample, priorSample);
            auto wristOrientationSample =
                EgocentricWristOrientation<Handedness>::TryCreate(sample, priorSample);
            if (handShapeSample.has_value() && wristOrientationSample.has_value())
            {
                return HandPose{ std::move(handShapeSample.value()), std::move(wristOrientationSample.value()) };
            }
            return{};
        }

        static NumberT Distance(const HandPose& a, const HandPose& b, gsl::span<const NumberT> tuning)
        {
            auto handShapeDistance = HandShape<Handedness>::Distance(
                a.m_handShapeSample,
                b.m_handShapeSample,
                TuningT::template getTuning<HandShape<Handedness>>(tuning));
            auto wristOrientationDistance = EgocentricWristOrientation<Handedness>::Distance(
                a.m_wristOrientationSample,
                b.m_wristOrientationSample,
                TuningT::template getTuning<EgocentricWristOrientation<Handedness>>(tuning));
            return std::max(handShapeDistance, wristOrientationDistance);
        }

        static HandPose Lerp(const HandPose& a, const HandPose& b, NumberT t)
        {
            auto handShape = HandShape<Handedness>::Lerp(a.m_handShapeSample, b.m_handShapeSample, t);
            auto wristOrientation = EgocentricWristOrientation<Handedness>::Lerp(a.m_wristOrientationSample, b.m_wristOrientationSample, t);
            return{ std::move(handShape), std::move(wristOrientation) };
        }

        template<typename ExamplesT>
        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(const ExamplesT& examples)
        {
            auto handShapeTuning = HandShape<Handedness>::CalculateTuning(examples);
            auto wristOrientationTuning = EgocentricWristOrientation<Handedness>::CalculateTuning(examples);
            return arrayConcat(handShapeTuning, wristOrientationTuning);
        }

        HandPose() = default;

    private:
        HandShape<Handedness> m_handShapeSample;
        EgocentricWristOrientation<Handedness> m_wristOrientationSample;

        HandPose(
            HandShape<Handedness> handShapeSample,
            EgocentricWristOrientation<Handedness> wristOrientationSample)
            : m_handShapeSample{ std::move(handShapeSample) }
            , m_wristOrientationSample{ std::move(wristOrientationSample) }
        {
        }
    };

    // TODO: Review theory of delta descriptors and assess how they work with delta sampling. Make sure creation isn't oscillating.
    template<Handedness Handedness>
    class WristRotation
    {
    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };
        static constexpr auto HANDEDNESS{ Handedness };

        static std::optional<WristRotation> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if (!sample.HmdPose.has_value())
            {
                return{};
            }

            if constexpr (Handedness == Handedness::LeftHanded)
            {
                if (sample.LeftWristPose.has_value() && priorSample.LeftWristPose.has_value())
                {
                    return WristRotation{ sample, priorSample };
                }
            }
            else if constexpr (Handedness == Handedness::RightHanded)
            {
                if (sample.RightWristPose.has_value() && priorSample.RightWristPose.has_value())
                {
                    return WristRotation{ sample, priorSample };
                }
            }
            return{};
        }

        static NumberT Distance(const WristRotation& a, const WristRotation& b, gsl::span<const NumberT> tuning)
        {
            auto distance = a.m_deltaOrientation.angularDistance(b.m_deltaOrientation);
            return normalizeDistance(distance, tuning[0]);
        }

        static WristRotation Lerp(const WristRotation& a, const WristRotation& b, NumberT t)
        {
            QuaternionT aQuat{ a.m_deltaOrientation };
            QuaternionT bQuat{ b.m_deltaOrientation };
            WristRotation result{};
            result.m_deltaOrientation = aQuat.slerp(t, bQuat);
            return result;
        }

        template<typename ExamplesT>
        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(const ExamplesT& examples)
        {
            if (examples.size() < 2)
            {
                return DEFAULT_TUNING;
            }

            auto sequences = createDescriptorSequencesFromExamples<WristRotation<Handedness>>(examples, DEFAULT_TUNING);
            auto extendedSequences = createDescriptorSequencesFromExamples<WristRotation<Handedness>>(examples, DEFAULT_TUNING, 1.);

            auto tuning = DEFAULT_TUNING;
            NumberT maxConnectionCost = 0;
            for (const auto& extendedSequence : extendedSequences)
            {
                for (const auto& sequence : sequences)
                {
                    auto distanceFunction = [](const auto& a, const auto& b) {
                        return a.m_deltaOrientation.angularDistance(b.m_deltaOrientation);
                    };
                    auto result = DynamicTimeWarping::Distance<const WristRotation<Handedness>>(extendedSequence, sequence, distanceFunction);
                    maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                }
            }
            // This descriptor is comparatively noisy and should mostly serve as a failsafe for other descriptors, so we tune 
            // to its identicality threshold rather than its midpoint.
            tuning[0] = std::max(maxConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
            return tuning;
        }

        WristRotation() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.1 };
        static inline constexpr NumberT IRRECONCILABILITY_THRESHOLD{ 0.2 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD, IRRECONCILABILITY_THRESHOLD) };
        trivial::Quaternion m_deltaOrientation{};

        WristRotation(const InputSample& sample, const InputSample& priorSample)
        {
            constexpr bool isLeftHanded = Handedness == Handedness::LeftHanded;
            auto wristPose = getWristPose<Handedness>(sample);
            auto priorWristPose = getWristPose<Handedness>(priorSample);

            AngleAxisT angleAxis{ priorWristPose.rotation().inverse() * wristPose.rotation() };
            m_deltaOrientation = QuaternionT{ angleAxis };
        }
    };

    // TODO: Review theory of delta descriptors and assess how they work with delta sampling. Make sure creation isn't oscillating.
    template<Handedness Handedness>
    class EgocentricWristTranslation
    {
    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };
        static constexpr auto HANDEDNESS{ Handedness };

        static std::optional<EgocentricWristTranslation> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if (!sample.HmdPose.has_value())
            {
                return {};
            }

            if constexpr (Handedness == Handedness::LeftHanded)
            {
                if (sample.LeftWristPose.has_value() && priorSample.LeftWristPose.has_value())
                {
                    return EgocentricWristTranslation{ sample, priorSample };
                }
            }
            else if constexpr (Handedness == Handedness::RightHanded)
            {
                if (sample.RightWristPose.has_value() && priorSample.RightWristPose.has_value())
                {
                    return EgocentricWristTranslation{ sample, priorSample };
                }
            }
            return{};
        }

        static NumberT Distance(
            const EgocentricWristTranslation& a,
            const EgocentricWristTranslation& b,
            gsl::span<const NumberT> tuning)
        {
            auto distance = a.m_egocentricTemporalPosition.distance(b.m_egocentricTemporalPosition);
            return normalizeDistance(distance, tuning[0]);
        }

        static EgocentricWristTranslation Lerp(
            const EgocentricWristTranslation& a,
            const EgocentricWristTranslation& b,
            NumberT t)
        {
            VectorT aVec{ a.m_egocentricTemporalPosition };
            VectorT bVec{ b.m_egocentricTemporalPosition };
            EgocentricWristTranslation result{};
            result.m_egocentricTemporalPosition = (static_cast<NumberT>(1) - t) * aVec + t * bVec;
            return result;
        }

        template<typename ExamplesT>
        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(const ExamplesT& examples)
        {
            if (examples.size() < 2)
            {
                return DEFAULT_TUNING;
            }

            auto sequences = createDescriptorSequencesFromExamples<EgocentricWristTranslation<Handedness>>(examples, DEFAULT_TUNING);
            auto extendedSequences = createDescriptorSequencesFromExamples<EgocentricWristTranslation<Handedness>>(examples, DEFAULT_TUNING, 1.);

            auto tuning = DEFAULT_TUNING;
            NumberT maxConnectionCost = 0;
            for (const auto& extendedSequence : extendedSequences)
            {
                for (const auto& sequence : sequences)
                {
                    auto distanceFunction = [](const auto& a, const auto& b) {
                        return a.m_egocentricTemporalPosition.distance(b.m_egocentricTemporalPosition);
                        };
                    auto result = DynamicTimeWarping::Match<const EgocentricWristTranslation<Handedness>>(extendedSequence, sequence, distanceFunction);
                    maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                }
            }
            NumberT midpoint = (IDENTICALITY_THRESHOLD + IRRECONCILABILITY_THRESHOLD) / 2;
            tuning[0] = std::max(maxConnectionCost / midpoint, DEFAULT_TUNING[0]);
            return tuning;
        }

        EgocentricWristTranslation() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.01 };
        static inline constexpr NumberT IRRECONCILABILITY_THRESHOLD{ 0.02 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD, IRRECONCILABILITY_THRESHOLD) };
        trivial::Point m_egocentricTemporalPosition{};

        EgocentricWristTranslation(const InputSample& sample, const InputSample& priorSample)
        {
            constexpr bool isLeftHanded = Handedness == Handedness::LeftHanded;
            auto& wristPose = isLeftHanded ? sample.LeftWristPose.value() : sample.RightWristPose.value();
            auto& priorWristPose =
                isLeftHanded ? priorSample.LeftWristPose.value() : priorSample.RightWristPose.value();
            auto& priorHmdPose = priorSample.HmdPose.value();

            auto ets = EgocentricTemporalSpace::getPose(priorWristPose.translation(), priorHmdPose);
            VectorT translation{ ets.inverse() * wristPose.translation() };
            m_egocentricTemporalPosition = translation;
        }
    };

    template<Handedness Handedness>
    class HandGesture
    {
    public:
        using TuningT = Tuning<HandPose<Handedness>, WristRotation<Handedness>, EgocentricWristTranslation<Handedness>>;
        static constexpr auto DEFAULT_TUNING{ TuningT::DEFAULT_TUNING };
        static constexpr auto HANDEDNESS{ Handedness };

        static std::optional<HandGesture> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample) {
            auto handPoseSample = HandPose<Handedness>::TryCreate(sample, priorSample);
            auto wristRotationSample = WristRotation<Handedness>::TryCreate(sample, priorSample);
            auto wristTranslationSample = EgocentricWristTranslation<Handedness>::TryCreate(sample, priorSample);
            if (handPoseSample.has_value() && wristRotationSample.has_value() && wristTranslationSample.has_value())
            {
                return HandGesture{
                    std::move(handPoseSample.value()), std::move(wristRotationSample.value()), std::move(wristTranslationSample.value()) };
            }
            return{};
        }

        static NumberT
            Distance(const HandGesture& a, const HandGesture& b, gsl::span<const NumberT> tuning)
        {
            auto handPoseDistance = HandPose<Handedness>::Distance(
                a.m_handPoseSample,
                b.m_handPoseSample,
                TuningT::template getTuning<HandPose<Handedness>>(tuning));
            auto wristRotationDistance = WristRotation<Handedness>::Distance(
                a.m_wristRotationSample,
                b.m_wristRotationSample,
                TuningT::template getTuning<WristRotation<Handedness>>(tuning));
            auto wristTranslationDistance = EgocentricWristTranslation<Handedness>::Distance(
                a.m_wristTranslationSample,
                b.m_wristTranslationSample,
                TuningT::template getTuning<EgocentricWristTranslation<Handedness>>(tuning));
            return std::max(handPoseDistance, std::max(wristRotationDistance, wristTranslationDistance));
        }

        static HandGesture Lerp(const HandGesture& a, const HandGesture& b, NumberT t)
        {
            auto handPose = HandPose<Handedness>::Lerp(a.m_handPoseSample, b.m_handPoseSample, t);
            auto wristRotation = WristRotation<Handedness>::Lerp(a.m_wristRotationSample, b.m_wristRotationSample, t);
            auto wristTranslation = EgocentricWristTranslation<Handedness>::Lerp(a.m_wristTranslationSample, b.m_wristTranslationSample, t);
            return{ std::move(handPose), std::move(wristRotation), std::move(wristTranslation) };
        }

        template<typename ExamplesT>
        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(const ExamplesT& examples)
        {
            auto handPoseTuning = HandPose<Handedness>::CalculateTuning(examples);
            auto wristRotationTuning = WristRotation<Handedness>::CalculateTuning(examples);
            auto wristTranslationTuning = EgocentricWristTranslation<Handedness>::CalculateTuning(examples);
            return arrayConcat(handPoseTuning, wristRotationTuning, wristTranslationTuning);
        }

        HandGesture() = default;

    private:
        HandPose<Handedness> m_handPoseSample;
        WristRotation<Handedness> m_wristRotationSample;
        EgocentricWristTranslation<Handedness> m_wristTranslationSample;

        HandGesture(
            HandPose<Handedness> handPoseSample,
            WristRotation<Handedness> wristRotationSample,
            EgocentricWristTranslation<Handedness> wristTranslationSample)
            : m_handPoseSample{ std::move(handPoseSample) }
            , m_wristRotationSample{ std::move(wristRotationSample) }
            , m_wristTranslationSample{ std::move(wristTranslationSample) }
        {
        }
    };

    class EgocentricRelativeWristPosition
    {
    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };
        static constexpr auto HANDEDNESS{ Handedness::TwoHanded };

        static std::optional<EgocentricRelativeWristPosition> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if (sample.HmdPose.has_value() &&
                sample.LeftWristPose.has_value() &&
                sample.RightWristPose.has_value())
            {
                return EgocentricRelativeWristPosition{ sample, priorSample };
            }
            return{};
        }

        static NumberT Distance(
            const EgocentricRelativeWristPosition& a,
            const EgocentricRelativeWristPosition& b,
            gsl::span<const NumberT> tuning)
        {
            auto distance = a.m_egocentricRelativeWristPosition.distance(b.m_egocentricRelativeWristPosition);
            return normalizeDistance(distance, tuning[0]);
        }

        static EgocentricRelativeWristPosition Lerp(
            const EgocentricRelativeWristPosition& a,
            const EgocentricRelativeWristPosition& b,
            NumberT t)
        {
            VectorT aVec{ a.m_egocentricRelativeWristPosition };
            VectorT bVec{ b.m_egocentricRelativeWristPosition };
            EgocentricRelativeWristPosition result{};
            result.m_egocentricRelativeWristPosition = (static_cast<NumberT>(1) - t) * aVec + t * bVec;
            return result;
        }

        template<typename ExamplesT>
        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(const ExamplesT& examples)
        {
            if (examples.size() < 2)
            {
                return DEFAULT_TUNING;
            }

            auto sequences = createDescriptorSequencesFromExamples<EgocentricRelativeWristPosition>(examples, DEFAULT_TUNING);
            auto extendedSequences = createDescriptorSequencesFromExamples<EgocentricRelativeWristPosition>(examples, DEFAULT_TUNING, 1.);

            auto tuning = DEFAULT_TUNING;
            NumberT maxConnectionCost = 0;
            for (const auto& extendedSequence : extendedSequences)
            {
                for (const auto& sequence : sequences)
                {
                    auto distanceFunction = [](const auto& a, const auto& b) {
                        return a.m_egocentricRelativeWristPosition.distance(b.m_egocentricRelativeWristPosition);
                        };
                    auto result = DynamicTimeWarping::Match<const EgocentricRelativeWristPosition>(extendedSequence, sequence, distanceFunction);
                    maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                }
            }
            NumberT midpoint = (IDENTICALITY_THRESHOLD + IRRECONCILABILITY_THRESHOLD) / 2;
            tuning[0] = std::max(maxConnectionCost / midpoint, DEFAULT_TUNING[0]);
            return tuning;
        }

        EgocentricRelativeWristPosition() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.04 };
        static inline constexpr NumberT IRRECONCILABILITY_THRESHOLD{ 0.2 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD, IRRECONCILABILITY_THRESHOLD) };
        trivial::Point m_egocentricRelativeWristPosition{};

        EgocentricRelativeWristPosition(const InputSample& sample, const InputSample&)
        {
            auto& leftWristPose = sample.LeftWristPose.value();
            auto& rightWristPose = sample.RightWristPose.value();
            auto& hmdPose = sample.HmdPose.value();

            auto ets = EgocentricTemporalSpace::getPose(leftWristPose.translation(), hmdPose);
            m_egocentricRelativeWristPosition = ets.inverse() * rightWristPose.translation();
        }
    };

    class TwoHandGesture {
    public:
        using TuningT = Tuning<
            HandGesture<Handedness::LeftHanded>,
            HandGesture<Handedness::RightHanded>,
            EgocentricRelativeWristPosition>;
        static constexpr auto DEFAULT_TUNING{ TuningT::DEFAULT_TUNING };
        static constexpr auto HANDEDNESS{ Handedness::TwoHanded };

        static std::optional<TwoHandGesture> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            auto leftGestureSample = HandGesture<Handedness::LeftHanded>::TryCreate(sample, priorSample);
            auto rightGestureSample = HandGesture<Handedness::RightHanded>::TryCreate(sample, priorSample);
            auto relativeSample = EgocentricRelativeWristPosition::TryCreate(sample, priorSample);
            if (leftGestureSample.has_value() &&
                rightGestureSample.has_value() &&
                relativeSample.has_value())
            {
                return TwoHandGesture{
                    std::move(leftGestureSample.value()),
                    std::move(rightGestureSample.value()),
                    std::move(relativeSample.value()) };
            }
            return{};
        }

        static NumberT Distance(const TwoHandGesture& a, const TwoHandGesture& b, gsl::span<const NumberT> tuning)
        {
            auto leftGestureDistance = HandGesture<Handedness::LeftHanded>::Distance(
                a.m_leftGestureSample,
                b.m_leftGestureSample,
                TuningT::template getTuning<HandGesture<Handedness::LeftHanded>>(tuning));
            auto rightGestureDistance = HandGesture<Handedness::RightHanded>::Distance(
                a.m_rightGestureSample,
                b.m_rightGestureSample,
                TuningT::template getTuning<HandGesture<Handedness::RightHanded>>(tuning));
            auto relativeWristPositionDistance = EgocentricRelativeWristPosition::Distance(
                a.m_relativeSample,
                b.m_relativeSample,
                TuningT::template getTuning<EgocentricRelativeWristPosition>(tuning));
            return std::max(leftGestureDistance, std::max(rightGestureDistance, relativeWristPositionDistance));
        }

        static TwoHandGesture Lerp(const TwoHandGesture& a, const TwoHandGesture& b, NumberT t)
        {
            auto leftGesture = HandGesture<Handedness::LeftHanded>::Lerp(a.m_leftGestureSample, b.m_leftGestureSample, t);
            auto rightGesture = HandGesture<Handedness::RightHanded>::Lerp(a.m_rightGestureSample, b.m_rightGestureSample, t);
            auto relativeWristPosition = EgocentricRelativeWristPosition::Lerp(a.m_relativeSample, b.m_relativeSample, t);
            return{ leftGesture, rightGesture, relativeWristPosition };
        }

        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(gsl::span<const action::Example> examples)
        {
            auto leftTuning = HandGesture<Handedness::LeftHanded>::CalculateTuning(examples);
            auto rightTuning = HandGesture<Handedness::RightHanded>::CalculateTuning(examples);
            auto relativeTuning = EgocentricRelativeWristPosition::CalculateTuning(examples);
            return arrayConcat(leftTuning, rightTuning, relativeTuning);
        }

        TwoHandGesture() = default;

    private:
        HandGesture<Handedness::LeftHanded> m_leftGestureSample;
        HandGesture<Handedness::RightHanded> m_rightGestureSample;
        EgocentricRelativeWristPosition m_relativeSample;

        TwoHandGesture(
            HandGesture<Handedness::LeftHanded> leftGestureSample,
            HandGesture<Handedness::RightHanded> rightGestureSample,
            EgocentricRelativeWristPosition relativeSample)
            : m_leftGestureSample{ std::move(leftGestureSample) }
            , m_rightGestureSample{ std::move(rightGestureSample) }
            , m_relativeSample{ std::move(relativeSample) }
        {
        }
    };

    template<typename DescriptorT>
    class TimestampedDescriptor
    {
    public:
        using TuningT = typename DescriptorT::TuningT;
        static constexpr auto DEFAULT_TUNING{ TuningT::DEFAULT_TUNING };
        static constexpr auto HANDEDNESS{ DescriptorT::HANDEDNESS };

        static std::optional<TimestampedDescriptor> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            auto underlyingDescriptor = DescriptorT::TryCreate(sample, priorSample);
            if (underlyingDescriptor.has_value())
            {
                return TimestampedDescriptor{ std::move(*underlyingDescriptor), sample.Timestamp };
            }
            return{};
        }

        static NumberT Distance(const TimestampedDescriptor& a, const TimestampedDescriptor& b, gsl::span<const NumberT> tuning)
        {
            return DescriptorT::Distance(a.m_underlyingDescriptor, b.m_underlyingDescriptor, tuning);
        }

        static TimestampedDescriptor Lerp(const TimestampedDescriptor& a, const TimestampedDescriptor& b, NumberT t)
        {
            auto underlyingDescriptor = DescriptorT::Lerp(a.m_underlyingDescriptor, b.m_underlyingDescriptor, t);
            auto timestamp = (static_cast<NumberT>(1) - t) * a.m_timestamp + t * b.m_timestamp;
            return{ std::move(underlyingDescriptor), timestamp };
        }

        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(gsl::span<const action::Example> examples)
        {
            return DescriptorT::CalculateTuning(examples);
        }

        TimestampedDescriptor() = default;

        const DescriptorT& getUnderlyingDescriptor() const
        {
            return m_underlyingDescriptor;
        }

        double getTimestamp()
        {
            return m_timestamp;
        }

    private:
        DescriptorT m_underlyingDescriptor{};
        double m_timestamp{};

        TimestampedDescriptor(DescriptorT underlyingDescriptor, double timestamp)
            : m_underlyingDescriptor{ std::move(underlyingDescriptor) }
            , m_timestamp{ timestamp }
        {}
    };
}
