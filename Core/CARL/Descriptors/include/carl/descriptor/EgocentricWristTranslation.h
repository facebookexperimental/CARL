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
    // TODO: Review theory of delta descriptors and assess how they work with delta sampling. Make sure creation isn't oscillating.
    template<Handedness Handedness>
    class EgocentricWristTranslation
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Egocentric Wrist Translation";

    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

        static NumberT AbsoluteDistance(
            const EgocentricWristTranslation& a,
            const EgocentricWristTranslation& b,
            gsl::span<const NumberT> tuning)
        {
            return Distance(a, a, b, b, tuning);
        }

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

            constexpr auto tryCreate = [](const InputSample& current, const InputSample& prior) { return TryCreate(current, prior); };
            auto sequences = createDescriptorSequencesFromExamples<EgocentricWristTranslation<Handedness>>(examples, DEFAULT_TUNING, tryCreate);
            auto extendedSequences = createDescriptorSequencesFromExamples<EgocentricWristTranslation<Handedness>>(examples, DEFAULT_TUNING, tryCreate, 1.);

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
                    auto result = DynamicTimeWarping::Match<const EgocentricWristTranslation<Handedness>>(extendedSequence, sequence, absoluteDistFn, deltaDistFn);
                    maxAverageConnectionCost = std::max<NumberT>(result.MaxConnectionCost / result.Connections, maxAverageConnectionCost);
                }
            }
            tuning[0] = std::max(maxAverageConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
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
            DynamicTimeWarping::Match<const EgocentricWristTranslation<Handedness>, decltype(absoluteDistFn), decltype(deltaDistFn), NumberT, true, decltype(rowsCallback)>(target, query, absoluteDistFn, deltaDistFn, 0, rowsCallback);
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
}
