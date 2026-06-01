/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/Recording.h>
#include <carl/descriptor/TimestampedDescriptor.h>
#include <carl/descriptor/DescriptorTraits.h>
#include <carl/descriptor/SequenceOperations.h>
#include <carl/DynamicTimeWarping.h>

#include <gsl/span>

#include <limits>
#include <vector>

namespace carl::action
{
    template<typename DescriptorT, typename CallableT>
    std::pair<double, double> autoTrimRecording(
        const Recording& recording,
        gsl::span<const std::vector<DescriptorT>> templates,
        gsl::span<const NumberT> tuning,
        const CallableT& tryCreate)
    {
        std::vector<descriptor::TimestampedDescriptor<DescriptorT>> timestampedSequence{};
        auto samples = recording.getSamples();
        auto mostRecentSample = samples.front();
        for (const auto& sample : samples)
        {
            descriptor::extendSequence<descriptor::TimestampedDescriptor<DescriptorT>>(
                sample, timestampedSequence, mostRecentSample, tuning, tryCreate);
        }

        std::vector<DescriptorT> sequence{};
        sequence.reserve(timestampedSequence.size());
        for (const auto& desc : timestampedSequence)
        {
            sequence.push_back(desc.getUnderlyingDescriptor());
        }

        auto absDist = [&tuning](const DescriptorT& a, const DescriptorT& b) {
            if constexpr (descriptor::DescriptorTraits<DescriptorT>::HAS_ABSOLUTE_DISTANCE)
            {
                return DescriptorT::AbsoluteDistance(a, b, tuning);
            }
            else
            {
                return NumberT{0};
            }
        };
        auto deltaDist = [&tuning](const DescriptorT& a, const DescriptorT& a0, const DescriptorT& b, const DescriptorT& b0) {
            if constexpr (descriptor::DescriptorTraits<DescriptorT>::HAS_DELTA_DISTANCE)
            {
                return DescriptorT::DeltaDistance(a, a0, b, b0, tuning);
            }
            else
            {
                return NumberT{0};
            }
        };

        auto bestMatchResult = DynamicTimeWarping::Match<const DescriptorT>(sequence, templates[0], absDist, deltaDist);
        for (size_t idx = 1; idx < templates.size(); ++idx)
        {
            auto matchResult = DynamicTimeWarping::Match<const DescriptorT>(sequence, templates[idx], absDist, deltaDist);
            if (matchResult.MatchCost < bestMatchResult.MatchCost)
            {
                bestMatchResult = matchResult;
            }
        }

        auto firstDescriptorT = timestampedSequence[bestMatchResult.ImageStartIdx].getTimestamp();
        auto lastDescriptorT = timestampedSequence[bestMatchResult.ImageStartIdx + bestMatchResult.ImageSize - 1].getTimestamp();
        constexpr auto epsilonT = std::numeric_limits<double>::epsilon();

        auto startT = std::numeric_limits<NumberT>::lowest() + 2 * epsilonT;
        auto endT = std::numeric_limits<NumberT>::max() - 2 * epsilonT;
        for (const auto& sample : samples)
        {
            if (sample.Timestamp < firstDescriptorT && sample.Timestamp > startT)
            {
                startT = sample.Timestamp;
            }

            if (sample.Timestamp > lastDescriptorT && sample.Timestamp < endT)
            {
                endT = sample.Timestamp;
            }
        }

        return { startT - epsilonT, endT + epsilonT };
    }
}
