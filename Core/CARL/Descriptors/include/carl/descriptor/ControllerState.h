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
    template<Handedness Handedness>
    class ControllerState
    {
        static constexpr std::array<const char*, static_cast<size_t>(InputSample::ControllerInput::COUNT)> ANALYSIS_DIMENSION_NAMES{
            "PRIMARY_CLICK",
            "SECONDARY_CLICK",
            "THUMBSTICK_X",
            "THUMBSTICK_Y",
            "THUMBSTICK_CLICK",
            "SQUEEZE_VALUE",
            "TRIGGER_VALUE",
        };

    public:
        static constexpr std::array<NumberT, ANALYSIS_DIMENSION_NAMES.size()> DEFAULT_TUNING{ 1., 1., 1., 1., 1., 1., 1. };

        static std::optional<ControllerState> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if constexpr (Handedness == Handedness::LeftHanded)
            {
                if (sample.LeftControllerInput.has_value())
                {
                    return ControllerState{ sample.LeftControllerInput.value() };
                }
            }
            else if constexpr (Handedness == Handedness::RightHanded)
            {
                if (sample.RightControllerInput.has_value())
                {
                    return ControllerState{ sample.RightControllerInput.value() };
                }
            }
            return {};
        }

        static NumberT Distance(
            const ControllerState& a,
            const ControllerState&,
            const ControllerState& b,
            const ControllerState&,
            gsl::span<const NumberT> tuning)
        {
            NumberT distance = 0;
            for (size_t idx = 0; idx < a.m_positions.size(); ++idx) {
                distance += InternalNormalizedDistance(idx, a, b, tuning);
            }
            return distance;
        }

        static ControllerState Lerp(const ControllerState& a, const ControllerState& b, NumberT t)
        {
            ControllerState result{};
            for (size_t idx = 0; idx < DEFAULT_TUNING.size(); ++idx)
            {
                NumberT aVal{ a.m_positions[idx] };
                NumberT bVal{ b.m_positions[idx] };
                result.m_positions[idx] = (static_cast<NumberT>(1) - t) * aVal + t * bVal;
            }
            return result;
        }

        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(gsl::span<const action::Example> examples)
        {
            if (examples.size() < 2)
            {
                return DEFAULT_TUNING;
            }

            constexpr auto tryCreate = [](const InputSample& current, const InputSample& prior) { return TryCreate(current, prior); };
            auto sequences = createDescriptorSequencesFromExamples<ControllerState<Handedness>>(examples, DEFAULT_TUNING, tryCreate);
            auto extendedSequences = createDescriptorSequencesFromExamples<ControllerState<Handedness>>(examples, DEFAULT_TUNING, tryCreate, 1.);

            auto tuning = DEFAULT_TUNING;
            for (size_t idx = 0; idx < tuning.size(); ++idx)
            {
                NumberT maxAverageConnectionCost = 0;
                for (const auto& extendedSequence : extendedSequences)
                {
                    for (const auto& sequence : sequences)
                    {
                        auto distanceFunction = [idx](const auto& a, const auto&, const auto& b, const auto&) {
                            return InternalRawDistance(idx, a, b);
                            };
                        auto result = DynamicTimeWarping::Match<const ControllerState<Handedness>>(extendedSequence, sequence, distanceFunction);
                        maxAverageConnectionCost = std::max<NumberT>(result.MaxConnectionCost / result.Connections, maxAverageConnectionCost);
                    }
                }
                tuning[idx] = std::max(maxAverageConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[idx]);
            }
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const ControllerState<Handedness>> target,
            gsl::span<const ControllerState<Handedness>> query,
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
                DynamicTimeWarping::Match<const ControllerState<Handedness>, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, 0, rowsCallback);
            }
            return results;
        }

        ControllerState() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.005 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
        std::array<NumberT, DEFAULT_TUNING.size()> m_positions{};

        ControllerState(const std::array<NumberT, DEFAULT_TUNING.size()>& state)
        {
            m_positions = state;
        }

        static NumberT InternalRawDistance(size_t dimension, const ControllerState& a, const ControllerState& b)
        {
            return static_cast<NumberT>(std::pow(a.m_positions[dimension] - b.m_positions[dimension], 2));
        }

        static NumberT InternalNormalizedDistance(size_t dimension, const ControllerState& a, const ControllerState& b, gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(dimension, a, b), tuning[dimension]);
        }
    };
}