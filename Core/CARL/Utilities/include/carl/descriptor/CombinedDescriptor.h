/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/descriptor/DescriptorTraits.h>
#include <carl/descriptor/DescriptorUtils.h>
#include <carl/descriptor/Trivial.h>
#include <carl/descriptor/Tuning.h>

namespace carl::descriptor
{
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

        static NumberT AbsoluteDistance(
            const CombinedDescriptor& a,
            const CombinedDescriptor& b,
            gsl::span<const NumberT> tuning)
        {
            return InternalAbsoluteDistance<0, Ts...>(a, b, tuning);
        }

        static NumberT DeltaDistance(
            const CombinedDescriptor& a,
            const CombinedDescriptor& a0,
            const CombinedDescriptor& b,
            const CombinedDescriptor& b0,
            gsl::span<const NumberT> tuning)
        {
            return InternalDeltaDistance<0, Ts...>(a, a0, b, b0, tuning);
        }

        static NumberT Distance(
            const CombinedDescriptor& a,
            const CombinedDescriptor& a0,
            const CombinedDescriptor& b,
            const CombinedDescriptor& b0,
            gsl::span<const NumberT> tuning)
        {
            return AbsoluteDistance(a, b, tuning) + DeltaDistance(a, a0, b, b0, tuning);
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
                auto absoluteDistFn = [tuning](const auto& a, const auto& b) {
                    return AbsoluteDistance(a, b, tuning);
                    };
                auto deltaDistFn = [tuning](const auto& a, const auto& a0, const auto& b, const auto& b0) {
                    return DeltaDistance(a, a0, b, b0, tuning);
                    };
                auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
                DynamicTimeWarping::Match<const CombinedDescriptor, decltype(absoluteDistFn), decltype(deltaDistFn), NumberT, true, decltype(rowsCallback)>(target, query, absoluteDistFn, deltaDistFn, 0, rowsCallback);
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
        static NumberT InternalAbsoluteDistance(
            const CombinedDescriptor& a,
            const CombinedDescriptor& b,
            gsl::span<const NumberT> tuning)
        {
            NumberT distance = 0;
            if constexpr (DescriptorTraits<T>::HAS_ABSOLUTE_DISTANCE)
            {
                const T& aDesc = a.m_underlyingDescriptors.template get<Idx>();
                const T& bDesc = b.m_underlyingDescriptors.template get<Idx>();
                distance = T::AbsoluteDistance(aDesc, bDesc, TuningT::template getTuning<T>(tuning));
            }
            if constexpr (sizeof...(RemainderT) > 0)
            {
                distance += InternalAbsoluteDistance<Idx + 1, RemainderT...>(a, b, tuning);
            }
            return distance;
        }

        template<size_t Idx, typename T, typename... RemainderT>
        static NumberT InternalDeltaDistance(
            const CombinedDescriptor& a,
            const CombinedDescriptor& a0,
            const CombinedDescriptor& b,
            const CombinedDescriptor& b0,
            gsl::span<const NumberT> tuning)
        {
            NumberT distance = 0;
            if constexpr (DescriptorTraits<T>::HAS_DELTA_DISTANCE)
            {
                const T& aDesc = a.m_underlyingDescriptors.template get<Idx>();
                const T& a0Desc = a0.m_underlyingDescriptors.template get<Idx>();
                const T& bDesc = b.m_underlyingDescriptors.template get<Idx>();
                const T& b0Desc = b0.m_underlyingDescriptors.template get<Idx>();
                distance = T::DeltaDistance(aDesc, a0Desc, bDesc, b0Desc, TuningT::template getTuning<T>(tuning));
            }
            if constexpr (sizeof...(RemainderT) > 0)
            {
                distance += InternalDeltaDistance<Idx + 1, RemainderT...>(a, a0, b, b0, tuning);
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

    // TODO: Revise how this is done as part of the work to break up Descriptor.h
    //#define ENABLE_TIMESTAMPED_ANALYSIS
#ifdef ENABLE_TIMESTAMPED_ANALYSIS
    template<typename... Ts>
    using CombinedDescriptorT = TimestampedDescriptor<CombinedDescriptor<Ts...>>;
#else
    template<typename... Ts>
    using CombinedDescriptorT = CombinedDescriptor<Ts...>;
#endif

}
