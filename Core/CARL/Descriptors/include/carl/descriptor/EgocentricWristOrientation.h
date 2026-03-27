/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/descriptor/DescriptorUtils.h>
#include <carl/descriptor/EgocentricTemporalSpace.h>
#include <carl/descriptor/SequenceOperations.h>
#include <carl/descriptor/Trivial.h>

namespace carl::descriptor
{
    template<Handedness Handedness>
    class EgocentricWristOrientation
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Egocentric Wrist Orientation";

    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

        static NumberT AbsoluteDistance(
            const EgocentricWristOrientation& a,
            const EgocentricWristOrientation& b,
            gsl::span<const NumberT> tuning)
        {
            return Distance(a, a, b, b, tuning);
        }

        static std::optional<EgocentricWristOrientation> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if (!sample.HmdPose.has_value())
            {
                return{};
            }

            auto optionalWristPose = getWristPose<Handedness>(sample);
            if constexpr (Handedness == Handedness::LeftHanded)
            {
                if (optionalWristPose.has_value() && priorSample.LeftWristPose.has_value())
                {
                    return EgocentricWristOrientation{ sample, optionalWristPose.value() };
                }
            }
            else if constexpr (Handedness == Handedness::RightHanded)
            {
                if (optionalWristPose.has_value() && priorSample.RightWristPose.has_value())
                {
                    return EgocentricWristOrientation{ sample, optionalWristPose.value() };
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


            constexpr auto tryCreate = [](const InputSample& current, const InputSample& prior) { return TryCreate(current, prior); };
            auto sequences = createDescriptorSequencesFromExamples<EgocentricWristOrientation<Handedness>>(examples, DEFAULT_TUNING, tryCreate);
            auto extendedSequences = createDescriptorSequencesFromExamples<EgocentricWristOrientation<Handedness>>(examples, DEFAULT_TUNING, tryCreate, 1.);

            auto tuning = DEFAULT_TUNING;
            NumberT maxAverageConnectionCost = 0;
            for (const auto& extendedSequence : extendedSequences)
            {
                for (const auto& sequence : sequences)
                {
                    auto absoluteDistFn = [](const auto& a, const auto& b) {
                        return InternalRawDistance(a, b);
                        };
                    auto deltaDistFn = [](const auto&, const auto&, const auto&, const auto&) -> NumberT {
                        return 0;
                        };
                    auto result = DynamicTimeWarping::Match<const EgocentricWristOrientation<Handedness>>(extendedSequence, sequence, absoluteDistFn, deltaDistFn);
                    maxAverageConnectionCost = std::max<NumberT>(result.MaxConnectionCost / result.Connections, maxAverageConnectionCost);
                }
            }
            tuning[0] = std::max(maxAverageConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
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
            auto absoluteDistFn = [tuning](const auto& a, const auto& b) {
                if constexpr (NormalizeDistance)
                {
                    return InternalNormalizedDistance(a, b, tuning);
                }
                else
                {
                    return InternalRawDistance(a, b);
                }
                };
            auto deltaDistFn = [](const auto&, const auto&, const auto&, const auto&) -> NumberT {
                return 0;
                };
            auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
            DynamicTimeWarping::Match<const EgocentricWristOrientation<Handedness>, decltype(absoluteDistFn), decltype(deltaDistFn), NumberT, true, decltype(rowsCallback)>(target, query, absoluteDistFn, deltaDistFn, 0, rowsCallback);
            return results;
        }

        EgocentricWristOrientation() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.2 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
        trivial::Quaternion m_egocentricTemporalOrientation{};

        EgocentricWristOrientation(const InputSample& sample, const TransformT& wristPose)
        {
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
}
