/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <gsl/span>

#include <vector>

namespace carl::DynamicTimeWarping
{
    namespace internal
    {
        template<typename ArgT, typename CallbackT, typename... Ts>
        void invokeVariadicCallback(ArgT&& arg, CallbackT& callback, ...)
        {
            callback(std::forward<ArgT>(arg));
        }
    }

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

    template<typename NumberT = double>
    struct MatchResult
    {
        NumberT MatchCost{};
        size_t Connections{};
        NumberT MaxConnectionCost{};
        size_t ImageStartIdx{};
        size_t ImageSize{};
    };

    template <typename VectorT, typename CallableT, typename NumberT = double, bool ReturnAllResults = false, typename... Ts>
    MatchResult<NumberT> Match(gsl::span<VectorT> target, gsl::span<VectorT> query, CallableT& distance, size_t minimumImageIdx = 0, Ts&... ts)
    {
        struct Entry
        {
            NumberT Cost{};
            size_t Connections{};
            NumberT MaxConnectionCost{};
            size_t StartIdx{};
        };

        constexpr auto analyzeRow = [](gsl::span<const Entry> row) {
            std::vector<MatchResult<NumberT>> analysis{};
            analysis.reserve(row.size() - 1);
            for (size_t idx = 1; idx < row.size(); ++idx)
            {
                const auto& match = row[idx];
                MatchResult<NumberT> result{
                    match.Cost,
                    match.Connections,
                    match.MaxConnectionCost,
                    match.StartIdx,
                    idx - match.StartIdx
                };
                analysis.push_back(result);
            }
            return analysis;
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

        if constexpr (ReturnAllResults)
        {
            internal::invokeVariadicCallback(analyzeRow(currentRow), ts...);
        }

        for (size_t j = 1; j < query.size(); ++j)
        {
            priorRow.swap(currentRow);

            for (size_t i = 0; i < target.size(); ++i)
            {
                const auto& ul = priorRow[i];
                const auto& u = priorRow[i + 1];
                const auto& l = currentRow[i];

                ulCost = distance(target[i], target[ul.StartIdx],  query[j], query[0]);
                uCost = distance(target[i], target[u.StartIdx], query[j], query[0]);
                lCost = distance(target[i], target[l.StartIdx], query[j], query[0]);

                auto ancestor = ul;
                cost = ulCost;
                if (uCost + u.Cost < cost + ancestor.Cost)
                {
                    ancestor = u;
                    cost = uCost;
                }
                if (lCost + l.Cost < cost + ancestor.Cost)
                {
                    ancestor = l;
                    cost = lCost;
                }
                currentRow[i + 1] = { cost + ancestor.Cost, ancestor.Connections + 1, std::max(cost, ancestor.MaxConnectionCost), ancestor.StartIdx };
            }

            if constexpr (ReturnAllResults)
            {
                internal::invokeVariadicCallback(analyzeRow(currentRow), ts...);
            }
        }

        NumberT matchCost = std::numeric_limits<NumberT>::max();
        size_t matchIdx = minimumImageIdx;
        for (size_t idx = minimumImageIdx; idx < currentRow.size(); ++idx)
        {
            if (currentRow[idx].Cost < matchCost)
            {
                matchIdx = idx;
                matchCost = currentRow[matchIdx].Cost;
            }
        }

        const auto& match = currentRow[matchIdx];
        return MatchResult<NumberT>{
            match.Cost,
            match.Connections,
            match.MaxConnectionCost,
            match.StartIdx,
            matchIdx - match.StartIdx,
        };
    }
}
