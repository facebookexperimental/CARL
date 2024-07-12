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
#include <utility>

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

        struct Transform
        {
            std::array<NumberT, 16> M{};

            Transform() = default;

            Transform(const TransformT& transform)
            {
                std::memcpy(M.data(), transform.matrix().data(), sizeof(M));
            }

            operator TransformT() const
            {
                auto m = M;
                return TransformT{ Eigen::Map<Eigen::Matrix4<NumberT>>{m.data()} };
            }

            Point operator *(const Point& p) const
            {
                return Point{
                    M[0] * p.X + M[4] * p.Y + M[8] * p.Z + M[12],
                    M[1] * p.X + M[5] * p.Y + M[9] * p.Z + M[13],
                    M[2] * p.X + M[6] * p.Y + M[10] * p.Z + M[14],
                };
            }
        };

        template<typename...>
        class Tuple;

        template<>
        class Tuple<> {};

        template<typename T, typename... Ts>
        class Tuple<T, Ts...>
        {
        public:
            Tuple() = default;

            Tuple(T&& t, Ts&&... ts)
                : m_t{ std::forward<T>(t) }, m_tail{ std::forward<Ts>(ts)... }
            {}

            template<size_t Idx, typename = std::enable_if_t<Idx != 0 && Idx <= sizeof...(Ts)>>
            auto& get()
            {
                return m_tail.template get<Idx - 1>();
            }

            template<size_t Idx, typename = std::enable_if_t<Idx != 0 && Idx <= sizeof...(Ts)>>
            const auto& get() const
            {
                return m_tail.template get<Idx - 1>();
            }

            template<size_t Idx, typename = std::enable_if_t<Idx == 0>>
            T& get()
            {
                return m_t;
            }

            template<size_t Idx, typename = std::enable_if_t<Idx == 0>>
            const T& get() const
            {
                return m_t;
            }

        private:
            T m_t;
            Tuple<Ts...> m_tail;
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

    constexpr NumberT NULL_TUNING{ 1000000 };
    constexpr auto createDistanceNormalizationFunction(NumberT identicalityThreshold, NumberT irreconcilabilityThreshold)
    {
        return [identicalityThreshold, irreconcilabilityThreshold](NumberT distance, NumberT tuning) {
            auto lowerBound = tuning * identicalityThreshold;
            auto upperBound = tuning * irreconcilabilityThreshold;
            return std::max<NumberT>(0, std::pow<NumberT>((distance - lowerBound) / (upperBound - lowerBound), 3));
        };
    }

    constexpr auto createDistanceNormalizationFunction(NumberT identicalityThreshold)
    {
        return createDistanceNormalizationFunction(identicalityThreshold, 2 * identicalityThreshold);
    }

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

                auto outerDistance = DescriptorT::Distance(*sampleDesc, sequence.back(), sequence.back(), sequence.back(), tuning);
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
                    auto distance = DescriptorT::Distance(intermediateDesc, sequence.back(), sequence.back(), sequence.back(), tuning);
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
        while (idx < samples.size() - 1 && samples[idx].Timestamp < startTimestamp)
        {
            ++idx;
        }
        while (idx < samples.size() - 1 && samples[idx].Timestamp < endTimestamp)
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

    using AnalysisT = std::tuple<std::string, NumberT, NumberT, std::vector<std::vector<DynamicTimeWarping::MatchResult<NumberT>>>>;

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
        static constexpr std::array<const char*, JOINTS.size()> ANALYSIS_DIMENSION_NAMES{
            "Thumb",
            "Index Finger",
            "Middle Finger",
            "Ring Finger",
            "Little Finger",
        };
        static constexpr size_t NORMALIZATION_JOINT{
            static_cast<size_t>(InputSample::Joint::IndexFingerBase) };

    public:
        static constexpr std::array<NumberT, JOINTS.size()> DEFAULT_TUNING{ 1., 1., 1., 1., 1. };

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

        static NumberT Distance(
            const HandShape& a,
            const HandShape&,
            const HandShape& b,
            const HandShape&,
            gsl::span<const NumberT> tuning)
        {
            NumberT distance = 0;
            for (size_t idx = 0; idx < a.m_positions.size(); ++idx) {
                distance += InternalNormalizedDistance(idx, a, b, tuning);
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
                        auto distanceFunction = [idx](const auto& a, const auto&, const auto& b, const auto&) {
                            return InternalRawDistance(idx, a, b);
                        };
                        auto result = DynamicTimeWarping::Match<const HandShape<Handedness>>(extendedSequence, sequence, distanceFunction);
                        maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                    }
                }
                tuning[idx] = std::max(maxConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[idx]);
            }
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const HandShape<Handedness>> target, 
            gsl::span<const HandShape<Handedness>> query, 
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            for (size_t idx = 0; idx < DEFAULT_TUNING.size(); ++idx)
            {
                results[idx] = { ANALYSIS_DIMENSION_NAMES[idx], IDENTICALITY_THRESHOLD, tuning[idx], {} };
                auto& rows = std::get<3>(results[idx]);
                auto distanceFunction = [idx, tuning](const auto& a, const auto&, const auto& b, const auto&) {
                    if constexpr (NormalizeDistance)
                    {
                        return InternalNormalizedDistance(idx, a, b, tuning);
                    }
                    else
                    {
                        return InternalRawDistance(idx, a, b);
                    }
                };
                auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
                DynamicTimeWarping::Match<const HandShape<Handedness>, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, rowsCallback);
            }
            return results;
        }

        HandShape() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.01 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
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

        static NumberT InternalRawDistance(size_t dimension, const HandShape& a, const HandShape& b)
        {
            return a.m_positions[dimension].distance(b.m_positions[dimension]);
        }

        static NumberT InternalNormalizedDistance(size_t dimension, const HandShape& a, const HandShape& b, gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(dimension, a, b), tuning[dimension]);
        }
    };

    template<Handedness Handedness>
    class EgocentricWristOrientation
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Egocentric Wrist Orientation";

    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

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
            const EgocentricWristOrientation&,
            const EgocentricWristOrientation& b,
            const EgocentricWristOrientation&,
            gsl::span<const NumberT> tuning)
        {
            return InternalNormalizedDistance(a, b, tuning);
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
                    auto distanceFunction = [](const auto& a, const auto&, const auto& b, const auto&) {
                        return InternalRawDistance(a, b);
                    };
                    auto result = DynamicTimeWarping::Match<const EgocentricWristOrientation<Handedness>>(extendedSequence, sequence, distanceFunction);
                    maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                }
            }
            tuning[0] = std::max(maxConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const EgocentricWristOrientation<Handedness>> target,
            gsl::span<const EgocentricWristOrientation<Handedness>> query,
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            results[0] = { ANALYSIS_DIMENSION_NAME, IDENTICALITY_THRESHOLD, tuning.front(), {} };
            auto& rows = std::get<3>(results[0]);
            auto distanceFunction = [tuning](const auto& a, const auto&, const auto& b, const auto&) {
                if constexpr (NormalizeDistance)
                {
                    return InternalNormalizedDistance(a, b, tuning);
                }
                else
                {
                    return InternalRawDistance(a, b);
                }
            };
            auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
            DynamicTimeWarping::Match<const EgocentricWristOrientation<Handedness>, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, rowsCallback);
            return results;
        }

        EgocentricWristOrientation() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.2 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
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

        static NumberT InternalRawDistance(const EgocentricWristOrientation& a, const EgocentricWristOrientation& b)
        {
            return a.m_egocentricTemporalOrientation.angularDistance(b.m_egocentricTemporalOrientation);
        }

        static NumberT InternalNormalizedDistance(const EgocentricWristOrientation& a, const EgocentricWristOrientation& b, gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(a, b), tuning.front());
        }
    };

    // TODO: Review theory of delta descriptors and assess how they work with delta sampling. Make sure creation isn't oscillating.
    template<Handedness Handedness>
    class WristRotation
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Wrist Rotation";

    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

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

        static NumberT Distance(
            const WristRotation& a,
            const WristRotation&,
            const WristRotation& b,
            const WristRotation&,
            gsl::span<const NumberT> tuning)
        {
            return InternalNormalizedDistance(a, b, tuning);
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
                    auto distanceFunction = [](const auto& a, const auto&, const auto& b, const auto&) {
                        return InternalRawDistance(a, b);
                    };
                    auto result = DynamicTimeWarping::Match<const WristRotation<Handedness>>(extendedSequence, sequence, distanceFunction);
                    maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                }
            }
            tuning[0] = std::max(maxConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const WristRotation<Handedness>> target,
            gsl::span<const WristRotation<Handedness>> query,
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            results[0] = { ANALYSIS_DIMENSION_NAME, IDENTICALITY_THRESHOLD, tuning[0], {} };
            auto& rows = std::get<3>(results[0]);
            auto distanceFunction = [tuning](const auto& a, const auto&, const auto& b, const auto&) {
                if constexpr (NormalizeDistance)
                {
                    return InternalNormalizedDistance(a, b, tuning);
                }
                else
                {
                    return InternalRawDistance(a, b);
                }
            };
            auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
            DynamicTimeWarping::Match<const WristRotation<Handedness>, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, rowsCallback);
            return results;
        }

        WristRotation() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.06 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
        trivial::Quaternion m_deltaOrientation{};

        WristRotation(const InputSample& sample, const InputSample& priorSample)
        {
            constexpr bool isLeftHanded = Handedness == Handedness::LeftHanded;
            auto wristPose = getWristPose<Handedness>(sample);
            auto priorWristPose = getWristPose<Handedness>(priorSample);

            AngleAxisT angleAxis{ priorWristPose.rotation().inverse() * wristPose.rotation() };
            m_deltaOrientation = QuaternionT{ angleAxis };
        }

        static NumberT InternalRawDistance(const WristRotation& a, const WristRotation& b)
        {
            return a.m_deltaOrientation.angularDistance(b.m_deltaOrientation);
        }

        static NumberT InternalNormalizedDistance(const WristRotation& a, const WristRotation& b, gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(a, b), tuning.front());
        }
    };

    // TODO: Review theory of delta descriptors and assess how they work with delta sampling. Make sure creation isn't oscillating.
    template<Handedness Handedness>
    class EgocentricWristTranslation
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Egocentric Wrist Translation";

    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

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
            const EgocentricWristTranslation&,
            const EgocentricWristTranslation& b,
            const EgocentricWristTranslation&,
            gsl::span<const NumberT> tuning)
        {
            return InternalNormalizedDistance(a, b, tuning);
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
                    auto distanceFunction = [](const auto& a, const auto&, const auto& b, const auto&) {
                        return InternalRawDistance(a, b);
                    };
                    auto result = DynamicTimeWarping::Match<const EgocentricWristTranslation<Handedness>>(extendedSequence, sequence, distanceFunction);
                    maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                }
            }
            tuning[0] = std::max(maxConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const EgocentricWristTranslation<Handedness>> target,
            gsl::span<const EgocentricWristTranslation<Handedness>> query,
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            results[0] = { ANALYSIS_DIMENSION_NAME, IDENTICALITY_THRESHOLD, tuning.front(), {} };
            auto& rows = std::get<3>(results[0]);
            auto distanceFunction = [tuning](const auto& a, const auto&, const auto& b, const auto&) {
                if constexpr (NormalizeDistance)
                {
                    return InternalNormalizedDistance(a, b, tuning);
                }
                else
                {
                    return InternalRawDistance(a, b);
                }
            };
            auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
            DynamicTimeWarping::Match<const EgocentricWristTranslation<Handedness>, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, rowsCallback);
            return results;
        }

        EgocentricWristTranslation() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.007 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
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

        static NumberT InternalRawDistance(const EgocentricWristTranslation<Handedness>& a, const EgocentricWristTranslation<Handedness>& b)
        {
            return a.m_egocentricTemporalPosition.distance(b.m_egocentricTemporalPosition);
        }

        static NumberT InternalNormalizedDistance(const EgocentricWristTranslation<Handedness>& a, const EgocentricWristTranslation<Handedness>& b, gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(a, b), tuning.front());
        }
    };

    template<Handedness Handedness>
    class EgocentricWristDisplacement
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Egocentric Wrist Displacement";

    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

        static std::optional<EgocentricWristDisplacement> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if (!sample.HmdPose.has_value())
            {
                return{};
            }

            if constexpr (Handedness == Handedness::LeftHanded)
            {
                if (sample.LeftWristPose.has_value())
                {
                    return EgocentricWristDisplacement{ sample, priorSample };
                }
            }
            else if constexpr (Handedness == Handedness::RightHanded)
            {
                if (sample.RightWristPose.has_value())
                {
                    return EgocentricWristDisplacement{ sample, priorSample };
                }
            }
            return{};
        }

        static NumberT Distance(
            const EgocentricWristDisplacement& a,
            const EgocentricWristDisplacement& a0,
            const EgocentricWristDisplacement& b,
            const EgocentricWristDisplacement& b0,
            gsl::span<const NumberT> tuning)
        {
            return InternalNormalizedDistance(a, a0, b, b0, tuning);
        }

        static EgocentricWristDisplacement Lerp(
            const EgocentricWristDisplacement& a,
            const EgocentricWristDisplacement& b,
            NumberT t)
        {
            EgocentricWristDisplacement result{};
            result.m_position = (1 - t) * VectorT{ a.m_position } + t * VectorT{ b.m_position };
            result.m_inverseEts = math::Lerp(a.m_inverseEts, b.m_inverseEts, t);
            return result;
        }

        template<typename ExamplesT>
        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(const ExamplesT& examples)
        {
            if (examples.size() < 2)
            {
                return DEFAULT_TUNING;
            }

            auto sequences = createDescriptorSequencesFromExamples<EgocentricWristDisplacement<Handedness>>(examples, DEFAULT_TUNING);
            auto extendedSequences = createDescriptorSequencesFromExamples<EgocentricWristDisplacement<Handedness>>(examples, DEFAULT_TUNING, 1.);

            auto tuning = DEFAULT_TUNING;
            NumberT maxConnectionCost = 0;
            for (const auto& extendedSequence : extendedSequences)
            {
                for (const auto& sequence : sequences)
                {
                    auto distanceFunction = [](const auto& a, const auto& a0, const auto& b, const auto& b0) {
                        return InternalRawDistance(a, a0, b, b0);
                    };
                    auto result = DynamicTimeWarping::Match<const EgocentricWristDisplacement<Handedness>>(extendedSequence, sequence, distanceFunction);
                    maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                }
            }
            tuning[0] = std::max(maxConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const EgocentricWristDisplacement<Handedness>> target,
            gsl::span<const EgocentricWristDisplacement<Handedness>> query,
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            results[0] = { ANALYSIS_DIMENSION_NAME, IDENTICALITY_THRESHOLD, tuning.front(), {} };
            auto& rows = std::get<3>(results[0]);
            auto distanceFunction = [tuning](const auto& a, const auto& a0, const auto& b, const auto& b0) {
                if constexpr (NormalizeDistance)
                {
                    return InternalNormalizedDistance(a, a0, b, b0, tuning);
                }
                else
                {
                    return InternalRawDistance(a, a0, b, b0);
                }
            };
            auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
            DynamicTimeWarping::Match<const EgocentricWristDisplacement<Handedness>, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, rowsCallback);
            return results;
        }

        EgocentricWristDisplacement() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.06 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
        trivial::Point m_position{};
        trivial::Transform m_inverseEts{};

        EgocentricWristDisplacement(const InputSample& sample, const InputSample&)
        {
            constexpr bool isLeftHanded = Handedness == Handedness::LeftHanded;
            auto& wristPose = isLeftHanded ? sample.LeftWristPose.value() : sample.RightWristPose.value();

            auto ets = EgocentricTemporalSpace::getPose(wristPose.translation(), sample.HmdPose.value());
            m_position = ets.translation();
            m_inverseEts = ets.inverse();
        }

        static NumberT InternalRawDistance(
            const EgocentricWristDisplacement& a,
            const EgocentricWristDisplacement& a0,
            const EgocentricWristDisplacement& b,
            const EgocentricWristDisplacement& b0)
        {
            auto aPos = a0.m_inverseEts * a.m_position;
            auto bPos = b0.m_inverseEts * b.m_position;
            return aPos.distance(bPos);
        }
        
        static NumberT InternalNormalizedDistance(
            const EgocentricWristDisplacement& a,
            const EgocentricWristDisplacement& a0,
            const EgocentricWristDisplacement& b,
            const EgocentricWristDisplacement& b0,
            gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(a, a0, b, b0), tuning.front());
        }
    };

    class EgocentricRelativeWristPosition
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Egocentric Relative Wrist Position";

    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

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
            const EgocentricRelativeWristPosition&,
            const EgocentricRelativeWristPosition& b,
            const EgocentricRelativeWristPosition&,
            gsl::span<const NumberT> tuning)
        {
            return InternalNormalizedDistance(a, b, tuning);
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
                    auto distanceFunction = [](const auto& a, const auto&, const auto& b, const auto&) {
                        return InternalRawDistance(a, b);
                    };
                    auto result = DynamicTimeWarping::Match<const EgocentricRelativeWristPosition>(extendedSequence, sequence, distanceFunction);
                    maxConnectionCost = std::max<NumberT>(result.MaxConnectionCost, maxConnectionCost);
                }
            }
            tuning[0] = std::max(maxConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const EgocentricRelativeWristPosition> target,
            gsl::span<const EgocentricRelativeWristPosition> query,
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            results[0] = { ANALYSIS_DIMENSION_NAME, IDENTICALITY_THRESHOLD, tuning.front(), {} };
            auto& rows = std::get<3>(results[0]);
            auto distanceFunction = [tuning](const auto& a, const auto&, const auto& b, const auto&) {
                if constexpr (NormalizeDistance)
                {
                    return InternalNormalizedDistance(a, b, tuning);
                }
                else
                {
                    return InternalRawDistance(a, b);
                }
            };
            auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
            DynamicTimeWarping::Match<const EgocentricRelativeWristPosition, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, rowsCallback);
            return results;
        }

        EgocentricRelativeWristPosition() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.06 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
        trivial::Point m_egocentricRelativeWristPosition{};

        EgocentricRelativeWristPosition(const InputSample& sample, const InputSample&)
        {
            auto& leftWristPose = sample.LeftWristPose.value();
            auto& rightWristPose = sample.RightWristPose.value();
            auto& hmdPose = sample.HmdPose.value();

            auto ets = EgocentricTemporalSpace::getPose(leftWristPose.translation(), hmdPose);
            m_egocentricRelativeWristPosition = ets.inverse() * rightWristPose.translation();
        }

        static NumberT InternalRawDistance(const EgocentricRelativeWristPosition& a, const EgocentricRelativeWristPosition& b)
        {
            return a.m_egocentricRelativeWristPosition.distance(b.m_egocentricRelativeWristPosition);
        }

        static NumberT InternalNormalizedDistance(const EgocentricRelativeWristPosition& a, const EgocentricRelativeWristPosition& b, gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(a, b), tuning.front());
        }
    };

    template<typename... Ts>
    class CombinedDescriptor
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Combined Descriptor";

    public:
        using TuningT = Tuning<Ts...>;
        static constexpr auto DEFAULT_TUNING{ TuningT::DEFAULT_TUNING };

        static std::optional<CombinedDescriptor> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            CombinedDescriptor descriptor{};
            if (descriptor.InternalTryCreate<0, Ts...>(sample, priorSample))
            {
                return descriptor;
            }
            else
            {
                return{};
            }
        }

        static NumberT Distance(
            const CombinedDescriptor& a, 
            const CombinedDescriptor& a0,
            const CombinedDescriptor& b, 
            const CombinedDescriptor& b0,
            gsl::span<const NumberT> tuning)
        {
            return InternalDistance<0, Ts...>(a, a0, b, b0, tuning);
        }

        static CombinedDescriptor Lerp(const CombinedDescriptor& a, const CombinedDescriptor& b, NumberT t)
        {
            CombinedDescriptor lerped{};
            lerped.InternalLerp<0, Ts...>(a, b, t);
            return lerped;
        }

        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(gsl::span<const action::Example> examples)
        {
            return InternalCalculateTuning<0, Ts...>(examples);
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const CombinedDescriptor> target,
            gsl::span<const CombinedDescriptor> query,
            gsl::span<const NumberT> tuning)
        {
            auto underlyingAnalysis = InternalAnalyze<0, NormalizeDistance, Ts...>(target, query, tuning);

            if constexpr (NormalizeDistance)
            {
                std::array<AnalysisT, 1> results{};
                results[0] = { ANALYSIS_DIMENSION_NAME, -1, -1, {} };
                auto& rows = std::get<3>(results[0]);
                auto distanceFunction = [tuning](const auto& a, const auto& a0, const auto& b, const auto& b0) {
                    return Distance(a, a0, b, b0, tuning);
                };
                auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
                DynamicTimeWarping::Match<const CombinedDescriptor, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, rowsCallback);
                return arrayConcat(underlyingAnalysis, results);
            }
            else
            {
                return underlyingAnalysis;
            }
        }

        CombinedDescriptor() = default;

    private:
        trivial::Tuple<Ts...> m_underlyingDescriptors;

        template<size_t Idx, typename T, typename... RemainderT>
        bool InternalTryCreate(const InputSample& sample, const InputSample& priorSample)
        {
            auto descriptor = T::TryCreate(sample, priorSample);
            if (descriptor.has_value())
            {
                m_underlyingDescriptors.template get<Idx>() = descriptor.value();
                if constexpr (sizeof...(RemainderT) == 0)
                {
                    return true;
                }
                else
                {
                    return InternalTryCreate<Idx + 1, RemainderT...>(sample, priorSample);
                }
            }
            else
            {
                return false;
            }
        }

        template<size_t Idx, typename T, typename... RemainderT>
        static NumberT InternalDistance(
            const CombinedDescriptor& a,
            const CombinedDescriptor& a0,
            const CombinedDescriptor& b,
            const CombinedDescriptor& b0,
            gsl::span<const NumberT> tuning)
        {
            const T& aDesc = a.m_underlyingDescriptors.template get<Idx>();
            const T& a0Desc = a0.m_underlyingDescriptors.template get<Idx>();
            const T& bDesc = b.m_underlyingDescriptors.template get<Idx>();
            const T& b0Desc = b0.m_underlyingDescriptors.template get<Idx>();
            auto distance = T::Distance(aDesc, a0Desc, bDesc, b0Desc, TuningT::template getTuning<T>(tuning));
            if constexpr (sizeof...(RemainderT) > 0)
            {
                distance += InternalDistance<Idx + 1, RemainderT...>(a, a0, b, b0, tuning);
            }
            return distance;
        }

        template<size_t Idx, typename T, typename... RemainderT>
        void InternalLerp(const CombinedDescriptor& a, const CombinedDescriptor& b, NumberT t)
        {
            const T& aDesc = a.m_underlyingDescriptors.template get<Idx>();
            const T& bDesc = b.m_underlyingDescriptors.template get<Idx>();
            m_underlyingDescriptors.template get<Idx>() = T::Lerp(aDesc, bDesc, t);
            if constexpr (sizeof...(RemainderT) > 0)
            {
                InternalLerp<Idx + 1, RemainderT...>(a, b, t);
            }
        }

        template<size_t Idx, typename T, typename... RemainderT>
        static auto InternalCalculateTuning(gsl::span<const action::Example> examples)
        {
            auto tuning = T::CalculateTuning(examples);
            if constexpr (sizeof...(RemainderT) > 0)
            {
                return arrayConcat(tuning, InternalCalculateTuning<Idx + 1, RemainderT...>(examples));
            }
            else
            {
                return tuning;
            }
        }

        template<size_t Idx, bool NormalizeDistance, typename T, typename... RemainderT>
        static auto InternalAnalyze(
            gsl::span<const CombinedDescriptor> target, 
            gsl::span<const CombinedDescriptor> query, 
            gsl::span<const NumberT> tuning)
        {
            std::vector<T> targetDescriptors{};
            targetDescriptors.reserve(target.size());
            for (const auto& combined : target)
            {
                targetDescriptors.push_back(combined.m_underlyingDescriptors.template get<Idx>());
            }

            std::vector<T> queryDescriptors{};
            queryDescriptors.reserve(query.size());
            for (const auto& combined : query)
            {
                queryDescriptors.push_back(combined.m_underlyingDescriptors.template get<Idx>());
            }

            auto descriptorTuning = TuningT::template getTuning<T>(tuning);
            
            auto results = T::template Analyze<NormalizeDistance>(targetDescriptors, queryDescriptors, descriptorTuning);
            if constexpr (sizeof...(RemainderT) > 0)
            {
                return arrayConcat(results, InternalAnalyze<Idx + 1, NormalizeDistance, RemainderT...>(target, query, tuning));
            }
            else
            {
                return results;
            }
        }
    };

    template<typename DescriptorT>
    class TimestampedDescriptor
    {
    public:
        using TuningT = typename DescriptorT::TuningT;
        static constexpr auto DEFAULT_TUNING{ TuningT::DEFAULT_TUNING };

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

        static NumberT Distance(
            const TimestampedDescriptor& a,
            const TimestampedDescriptor& a0,
            const TimestampedDescriptor& b,
            const TimestampedDescriptor& b0,
            gsl::span<const NumberT> tuning)
        {
            return DescriptorT::Distance(a.m_underlyingDescriptor, a0.m_underlyingDescriptor, b.m_underlyingDescriptor, b0.m_underlyingDescriptor, tuning);
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

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const TimestampedDescriptor> target, 
            gsl::span<const TimestampedDescriptor> query, 
            gsl::span<const NumberT> tuning)
        {
            std::vector<DescriptorT> targetDescriptors{};
            targetDescriptors.reserve(target.size());
            for (const auto& descriptor : target)
            {
                targetDescriptors.push_back(descriptor.m_underlyingDescriptor);
            }

            std::vector<DescriptorT> queryDescriptors{};
            queryDescriptors.reserve(query.size());
            for (const auto& descriptor : query)
            {
                queryDescriptors.push_back(descriptor.m_underlyingDescriptor);
            }

            return DescriptorT::Analyze<NormalizeDistance>(targetDescriptors, queryDescriptors, tuning);
        }

        TimestampedDescriptor() = default;

        const DescriptorT& getUnderlyingDescriptor() const
        {
            return m_underlyingDescriptor;
        }

        double getTimestamp() const
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

// TODO: Revise how this is done as part of the work to break up Descriptor.h
#define ENABLE_TIMESTAMPED_ANALYSIS
#ifdef ENABLE_TIMESTAMPED_ANALYSIS
    template<typename... Ts>
    using CombinedDescriptorT = TimestampedDescriptor<CombinedDescriptor<Ts...>>;
#else
    template<typename... Ts>
    using CombinedDescriptorT = CombinedDescriptor<Ts...>;
#endif

    template<Handedness Handedness>
    using HandPose = CombinedDescriptorT<HandShape<Handedness>, EgocentricWristOrientation<Handedness>>;

    template<Handedness Handedness>
    using HandGesture = CombinedDescriptorT<HandPose<Handedness>, WristRotation<Handedness>, EgocentricWristTranslation<Handedness>, EgocentricWristDisplacement<Handedness>>;

    using TwoHandGesture = CombinedDescriptorT<HandGesture<Handedness::LeftHanded>, HandGesture<Handedness::RightHanded>, EgocentricRelativeWristPosition>;

    template<Handedness Handedness>
    using ControllerGesture = CombinedDescriptorT<EgocentricWristOrientation<Handedness>, WristRotation<Handedness>, EgocentricWristTranslation<Handedness>, EgocentricWristDisplacement<Handedness>>;

    using TwoControllerGesture = CombinedDescriptorT<ControllerGesture<Handedness::LeftHanded>, ControllerGesture<Handedness::RightHanded>, EgocentricRelativeWristPosition>;
}
