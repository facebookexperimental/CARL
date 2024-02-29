/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <gsl/span>

#include <array>
#include <vector>

namespace carl::DynamicTimeWarping
{
    template <typename VectorT, typename CallableT, typename NumberT = double>
    NumberT Distance(gsl::span<VectorT> a, gsl::span<VectorT> b, CallableT& distance)
    {
        thread_local std::vector<NumberT> priorRow{};
        thread_local std::vector<NumberT> currentRow{};
        priorRow.resize(a.size() + 1);
        currentRow.resize(priorRow.size());

        currentRow[0] = static_cast<NumberT>(0);
        for (size_t idx = 1; idx < currentRow.size(); ++idx)
        {
            currentRow[idx] = std::numeric_limits<NumberT>::max();
        }

        for (size_t j = 0; j < b.size(); ++j)
        {
            currentRow.swap(priorRow);

            currentRow[0] = std::numeric_limits<NumberT>::max();
            for (size_t i = 0; i < a.size(); ++i)
            {
                NumberT cost = distance(a[i], b[j]);
                currentRow[i + 1] = cost + std::min(priorRow[i], std::min(priorRow[i + 1], currentRow[i]));
            }
        }

        return currentRow.back();
    }

    // To allow fudge factors on an unsegmented sequence, reverse time and allow fudge in the answer.
    // When you compare the sequences backwards, you effectively ask, "Did this template end now?"
    // Allow the query sequence to be longer than the template sequence, then find the minimimum
    // cost along the bottom row near the length of the template sequence. This will be the "injective"
    // cost of matching the template against the "most recent" part of the input sequence. Note that 
    // the "minimum distance along the bottom row" operation must be assessed as a NORMALIZED distance
    // (the distance divided by the number of matches, which equals the larger of the matched sequences) 
    // or else the algorithm will simply find the shortest image, not the best, since every additional 
    // element introduces a nonnegative absolute distance.
    template <typename VectorT, typename CallableT, typename NumberT = double, bool ReverseTime = true>
    std::tuple<NumberT, size_t> InjectiveDistanceAndImageSize(
        gsl::span<VectorT> longer,
        gsl::span<VectorT> shorter,
        CallableT& distance,
        NumberT minimumImageRatio = 0)
    {
        auto& a = longer;
        auto& b = shorter;

        thread_local std::vector<NumberT> priorRow{};
        thread_local std::vector<NumberT> currentRow{};
        priorRow.resize(a.size() + 1);
        currentRow.resize(priorRow.size());

        currentRow[0] = static_cast<NumberT>(0);
        for (size_t idx = 1; idx < priorRow.size(); ++idx)
        {
            currentRow[idx] = std::numeric_limits<NumberT>::max();
        }

        for (size_t j = 0; j < b.size(); ++j)
        {
            priorRow.swap(currentRow);

            currentRow[0] = std::numeric_limits<NumberT>::max();
            for (size_t i = 0; i < a.size(); ++i)
            {
                NumberT cost;
                if constexpr (ReverseTime)
                {
                    cost = distance(a[a.size() - i - 1], b[b.size() - j - 1]);
                }
                else
                {
                    cost = distance(a[i], b[j]);
                }
                currentRow[i + 1] = cost + std::min(priorRow[i], std::min(priorRow[i + 1], currentRow[i]));
            }
        }

        size_t disallowedImageLength = static_cast<size_t>(std::floor(minimumImageRatio * shorter.size()));
        NumberT minimum = std::numeric_limits<NumberT>::max();
        size_t minimumIdx = disallowedImageLength;
        for (size_t idx = disallowedImageLength; idx < currentRow.size(); ++idx)
        {
            size_t connectionsCount = std::max(idx, b.size());
            if (currentRow[idx] / connectionsCount < minimum)
            {
                minimumIdx = idx;
                minimum = currentRow[minimumIdx] / connectionsCount;
            }
        }
        return { minimum, minimumIdx };
    }

    template <typename VectorT, typename CallableT, typename NumberT = double, bool ReverseTime = true>
    std::tuple<NumberT, size_t> AdaptiveStartInjectiveDistanceAndImageSize(
        gsl::span<VectorT> longer,
        gsl::span<VectorT> shorter,
        CallableT& distance,
        NumberT minimumImageRatio = 0)
    {
        size_t maxAdaptiveStartSize = std::min(longer.size() - shorter.size(), shorter.size() / 2);
        size_t idx = 0;
        NumberT score = distance(longer[idx], shorter[0]);
        for (; idx < maxAdaptiveStartSize; ++idx)
        {
            NumberT nextScore = distance(longer[idx + 1], shorter[0]);
            if (nextScore > score)
            {
                break;
            }
            score = nextScore;
        }
        gsl::span<VectorT> newLonger = gsl::make_span<VectorT>(&longer[idx], longer.size() - idx);
        return InjectiveDistanceAndImageSize(newLonger, shorter, distance, minimumImageRatio);
    }
}
