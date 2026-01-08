/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "DescriptorUtils.h"

namespace carl::descriptor
{
    template<typename DescriptorT>
    class TimestampedDescriptor
    {
    public:
        static constexpr auto DEFAULT_TUNING{ DescriptorT::DEFAULT_TUNING };

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
        {
        }
    };
}
