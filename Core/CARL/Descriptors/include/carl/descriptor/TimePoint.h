/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/descriptor/DescriptorUtils.h>
#include <carl/descriptor/SequenceOperations.h>

namespace carl::descriptor
{
    class TimePoint
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "TimePoint";

    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

        static NumberT DeltaDistance(
            const TimePoint& a,
            const TimePoint& a0,
            const TimePoint& b,
            const TimePoint& b0,
            gsl::span<const NumberT> tuning)
        {
            return Distance(a, a0, b, b0, tuning);
        }

        static std::optional<TimePoint> TryCreate(
            const InputSample& sample,
            const InputSample&)
        {
            return TimePoint{ sample };
        }

        static NumberT Distance(
            const TimePoint& a,
            const TimePoint& a0,
            const TimePoint& b,
            const TimePoint& b0,
            gsl::span<const NumberT> tuning)
        {
            return InternalNormalizedDistance(a, a0, b, b0, tuning);
        }

        static TimePoint Lerp(const TimePoint& a, const TimePoint& b, NumberT t)
        {
            TimePoint timepoint{};
            timepoint.m_timestamp = (NumberT{ 1 } - t) * a.m_timestamp + t * b.m_timestamp;
            return timepoint;
        }

        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(gsl::span<const action::Example> examples)
        {
            if (examples.size() < 2)
            {
                return DEFAULT_TUNING;
            }

            constexpr auto tryCreate = [](const InputSample& current, const InputSample& prior) { return TryCreate(current, prior); };
            auto sequences = createDescriptorSequencesFromExamples<TimePoint>(examples, DEFAULT_TUNING, tryCreate);
            auto extendedSequences = createDescriptorSequencesFromExamples<TimePoint>(examples, DEFAULT_TUNING, tryCreate, 1.);

            auto tuning = DEFAULT_TUNING;
            for (size_t idx = 0; idx < tuning.size(); ++idx)
            {
                NumberT maxAverageConnectionCost = 0;
                for (const auto& extendedSequence : extendedSequences)
                {
                    for (const auto& sequence : sequences)
                    {
                        auto absoluteDistFn = [](const auto&, const auto&) -> NumberT {
                            return 0;
                            };
                        auto deltaDistFn = [](const auto& a, const auto& a0, const auto& b, const auto& b0) {
                            return InternalRawDistance(a, a0, b, b0);
                            };
                        auto result = DynamicTimeWarping::Match<const TimePoint>(extendedSequence, sequence, absoluteDistFn, deltaDistFn);
                        maxAverageConnectionCost = std::max<NumberT>(result.MaxConnectionCost / result.Connections, maxAverageConnectionCost);
                    }
                }
                tuning[idx] = std::max(maxAverageConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[idx]);
            }
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const TimePoint> target,
            gsl::span<const TimePoint> query,
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            results[0] = { ANALYSIS_DIMENSION_NAME, IDENTICALITY_THRESHOLD, tuning.front(), {} };
            auto& rows = std::get<3>(results[0]);
            auto absoluteDistFn = [](const auto&, const auto&) -> NumberT {
                return 0;
                };
            auto deltaDistFn = [tuning](const auto& a, const auto& a0, const auto& b, const auto& b0) {
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
            DynamicTimeWarping::Match<const TimePoint, decltype(absoluteDistFn), decltype(deltaDistFn), NumberT, true, decltype(rowsCallback)>(target, query, absoluteDistFn, deltaDistFn, 0, rowsCallback);
            return results;
        }

        TimePoint() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.1 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
        NumberT m_timestamp{};

        TimePoint(const InputSample& sample)
            : m_timestamp{ static_cast<NumberT>(sample.Timestamp) }
        {
        }

        static NumberT InternalRawDistance(const TimePoint& a, const TimePoint& a0, const TimePoint& b, const TimePoint& b0)
        {
            const auto deltaA = a.m_timestamp - a0.m_timestamp;
            const auto deltaB = b.m_timestamp - b0.m_timestamp;
            return std::abs(deltaA - deltaB);
        }

        static NumberT InternalNormalizedDistance(const TimePoint& a, const TimePoint& a0, const TimePoint& b, const TimePoint& b0, gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(a, a0, b, b0), tuning.front());
        }
    };
}
