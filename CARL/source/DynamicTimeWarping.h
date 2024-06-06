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
    auto Distance(gsl::span<VectorT> a, gsl::span<VectorT> b, CallableT& distance)
    {
        struct Entry
        {
            NumberT Cost{};
            size_t Connections{};
            NumberT MaxConnectionCost{};
        };

        thread_local std::vector<Entry> priorRow{};
        thread_local std::vector<Entry> currentRow{};
        priorRow.resize(a.size() + 1);
        currentRow.resize(priorRow.size());

        currentRow[0] = { static_cast<NumberT>(0), 0 };
        for (size_t idx = 1; idx < currentRow.size(); ++idx)
        {
            currentRow[idx] = { std::numeric_limits<NumberT>::max(), 0, std::numeric_limits<NumberT>::max() };
        }

        for (size_t j = 0; j < b.size(); ++j)
        {
            currentRow.swap(priorRow);

            currentRow[0] = { std::numeric_limits<NumberT>::max(), 0, std::numeric_limits<NumberT>::max() };
            for (size_t i = 0; i < a.size(); ++i)
            {
                NumberT cost = distance(a[i], b[j]);
                auto ancestor = priorRow[i];
                if (priorRow[i + 1].Cost < ancestor.Cost)
                {
                    ancestor = priorRow[i + 1];
                }
                if (currentRow[i].Cost < ancestor.Cost)
                {
                    ancestor = currentRow[i];
                }
                currentRow[i + 1] = { cost + ancestor.Cost, ancestor.Connections + 1, std::max(cost, ancestor.MaxConnectionCost) };
            }
        }

        return currentRow.back();
    }

    // To allow fudge factors on an unsegmented sequence, reverse time and allow fudge in the answer.
    // When you compare the sequences backwards, you effectively ask, "Did this template end now?"
    // Allow the query sequence to be longer than the template sequence, then find the minimimum
    // cost along the bottom row near the length of the template sequence. This will be the "injective"
    // cost of matching the template against the "most recent" part of the input sequence.
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
                NumberT cost{};
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
            if (currentRow[idx] < minimum)
            {
                minimumIdx = idx;
                minimum = currentRow[minimumIdx];
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

    template <typename VectorT, typename CallableT, typename NumberT = double>
    auto Match(
        gsl::span<VectorT> target,
        gsl::span<VectorT> query,
        CallableT& distance)
    {
        struct Entry
        {
            NumberT Cost{};
            size_t Connections{};
            NumberT MaxConnectionCost{};
            size_t StartIdx{};
        };

        thread_local std::vector<Entry> priorRow{};
        thread_local std::vector<Entry> currentRow{};
        priorRow.resize(target.size() + 1);
        currentRow.resize(priorRow.size());

        NumberT ulCost{};
        NumberT uCost{};
        NumberT lCost{};

        NumberT cost{};

        constexpr Entry sentinel{ std::numeric_limits<NumberT>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<NumberT>::max(), 0 };
        priorRow[0] = sentinel;
        currentRow[0] = sentinel;

        for (size_t i = 0; i < target.size(); ++i)
        {
            cost = distance(target[i], target[i], query[0], query[0]);
            currentRow[i + 1] = { cost, 1, cost, i };
        }

        for (size_t j = 1; j < query.size(); ++j)
        {
            priorRow.swap(currentRow);

            for (size_t i = 0; i < target.size(); ++i)
            {
                ulCost = priorRow[i].Cost;
                ulCost += distance(target[i], target[priorRow[i].StartIdx],  query[j], query[0]);

                uCost = priorRow[i + 1].Cost;
                uCost += distance(target[i], target[priorRow[i + 1].StartIdx], query[j], query[0]);

                lCost = currentRow[i].Cost;
                lCost += distance(target[i], target[currentRow[i].StartIdx], query[j], query[0]);

                auto ancestor = priorRow[i];
                cost = ulCost;
                if (uCost < cost)
                {
                    ancestor = priorRow[i + 1];
                    cost = uCost;
                }
                if (lCost < cost)
                {
                    ancestor = currentRow[i];
                    cost = lCost;
                }
                currentRow[i + 1] = { cost, ancestor.Connections + 1, std::max(cost, ancestor.MaxConnectionCost), ancestor.StartIdx };
            }
        }

        NumberT matchCost = std::numeric_limits<NumberT>::max();
        size_t matchIdx = 0;
        for (size_t idx = 0; idx < currentRow.size(); ++idx)
        {
            if (currentRow[idx].Cost < matchCost)
            {
                matchIdx = idx;
                matchCost = currentRow[matchIdx].Cost;
            }
        }

        auto& match = currentRow[matchIdx];
        struct
        {
            NumberT MatchCost{};
            size_t Connections{};
            NumberT MaxConnectionCost{};
            size_t ImageStartIdx{};
            size_t ImageSize{};
        } result{
            match.Cost,
            match.Connections,
            match.MaxConnectionCost,
            match.StartIdx,
            matchIdx - match.StartIdx,
        };
        return result;
    }
}
