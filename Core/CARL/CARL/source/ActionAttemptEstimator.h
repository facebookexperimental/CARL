/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/ActionAttemptEstimation.h>
#include <carl/DynamicTimeWarping.h>
#include <carl/RingBuffer.h>
#include <carl/descriptor/DescriptorTraits.h>

#include "SessionImpl.h"

#include <gsl/span>

#include <optional>
#include <vector>

namespace carl::action
{
    template<typename DescriptorT>
    class ActionAttemptEstimator
    {
    public:
        explicit ActionAttemptEstimator(ActionAttemptEstimationSettings settings, gsl::span<const NumberT> tuning)
            : m_settings{ std::move(settings) }
            , m_distanceHistory{ 5 }
        {
            std::copy(tuning.begin(), tuning.end(), m_tuning.begin());

            if (m_settings.precomputedSignificanceThreshold.has_value())
            {
                m_ambientEMA = *m_settings.precomputedSignificanceThreshold / m_settings.significanceMultiple;
                m_warmedUp = true;
            }
        }

        void update(
            double timestamp,
            const DynamicTimeWarping::MatchResult<NumberT>& matchResult,
            size_t bestTemplateIdx,
            gsl::span<const DescriptorT> sequence,
            gsl::span<const double> timestamps,
            Session::Impl& sessionImpl,
            NumberT score)
        {
            NumberT distance = matchResult.Connections > 0
                ? matchResult.MatchCost / matchResult.Connections
                : std::numeric_limits<NumberT>::max();

            m_distanceHistory.push_back(distance);
            updateAmbientBaseline(distance);

            if (score > 0.7)
            {
                discardOverlappingPending(matchResult.ImageStartIdx, matchResult.ImageSize);
                return;
            }

            if (!m_warmedUp)
            {
                return;
            }

            if (isSignificantMinimum())
            {
                NumberT minimumDistance = m_distanceHistory[m_distanceHistory.size() - 2];

                if (m_pendingSnapshot.has_value() && overlaps(m_pendingSnapshot->imageStartIdx, m_pendingSnapshot->imageSize, matchResult.ImageStartIdx, matchResult.ImageSize))
                {
                    if (minimumDistance < m_pendingSnapshot->distance)
                    {
                        captureSnapshot(timestamp, minimumDistance, bestTemplateIdx, matchResult, sequence, timestamps, sessionImpl);
                    }
                }
                else
                {
                    commitPendingSnapshot();
                    captureSnapshot(timestamp, minimumDistance, bestTemplateIdx, matchResult, sequence, timestamps, sessionImpl);
                }
            }
        }

        void expireSnapshots(double currentTimestamp)
        {
            auto it = std::remove_if(m_snapshots.begin(), m_snapshots.end(),
                [&](const Snapshot& s) { return (currentTimestamp - s.timestamp) > m_settings.snapshotExpirySeconds; });
            m_snapshots.erase(it, m_snapshots.end());

            if (m_pendingSnapshot.has_value() && (currentTimestamp - m_pendingSnapshot->timestamp) > m_settings.snapshotExpirySeconds)
            {
                m_pendingSnapshot.reset();
            }
        }

        std::optional<AttemptCluster> checkForCluster()
        {
            commitPendingSnapshot();

            if (m_snapshots.size() < m_settings.clusterMinCount)
            {
                return {};
            }

            double timeSpan = m_snapshots.back().timestamp - m_snapshots.front().timestamp;
            if (timeSpan > m_settings.clusterWindowSeconds)
            {
                return {};
            }

            if (!snapshotsCluster())
            {
                return {};
            }

            AttemptCluster cluster{};
            cluster.firstTimestamp = m_snapshots.front().timestamp;
            cluster.lastTimestamp = m_snapshots.back().timestamp;
            cluster.count = m_snapshots.size();
            for (auto& snapshot : m_snapshots)
            {
                cluster.snapshots.push_back(std::move(snapshot.inputSamples));
            }
            m_snapshots.clear();
            return cluster;
        }

    private:
        struct Snapshot
        {
            std::vector<DescriptorT> descriptors{};
            std::vector<InputSample> inputSamples{};
            double timestamp{};
            NumberT distance{};
            size_t imageStartIdx{};
            size_t imageSize{};
        };

        ActionAttemptEstimationSettings m_settings{};
        std::array<NumberT, DescriptorT::DEFAULT_TUNING.size()> m_tuning{};
        NumberT m_ambientEMA{};
        size_t m_updateCount{};
        bool m_warmedUp{};
        RingBuffer<NumberT> m_distanceHistory;
        std::optional<Snapshot> m_pendingSnapshot{};
        std::vector<Snapshot> m_snapshots{};

        static constexpr size_t WARMUP_FRAMES{ 60 };

        void updateAmbientBaseline(NumberT distance)
        {
            if (m_updateCount == 0)
            {
                m_ambientEMA = distance;
            }
            else
            {
                m_ambientEMA = m_settings.ambientAlpha * distance + (1.0 - m_settings.ambientAlpha) * m_ambientEMA;
            }
            ++m_updateCount;

            if (!m_warmedUp && m_updateCount >= WARMUP_FRAMES)
            {
                m_warmedUp = true;
            }
        }

        bool isSignificantMinimum() const
        {
            if (m_distanceHistory.size() < 3)
            {
                return false;
            }

            size_t sz = m_distanceHistory.size();
            NumberT prev = m_distanceHistory[sz - 3];
            NumberT candidate = m_distanceHistory[sz - 2];
            NumberT curr = m_distanceHistory[sz - 1];

            return candidate < prev && candidate < curr
                && candidate < m_ambientEMA * m_settings.significanceMultiple;
        }

        static bool overlaps(size_t startA, size_t sizeA, size_t startB, size_t sizeB)
        {
            size_t endA = startA + sizeA;
            size_t endB = startB + sizeB;
            return startA < endB && startB < endA;
        }

        void discardOverlappingPending(size_t imageStartIdx, size_t imageSize)
        {
            if (m_pendingSnapshot.has_value() && overlaps(m_pendingSnapshot->imageStartIdx, m_pendingSnapshot->imageSize, imageStartIdx, imageSize))
            {
                m_pendingSnapshot.reset();
            }
        }

        void commitPendingSnapshot()
        {
            if (m_pendingSnapshot.has_value())
            {
                m_snapshots.push_back(std::move(*m_pendingSnapshot));
                m_pendingSnapshot.reset();
            }
        }

        void captureSnapshot(
            double timestamp,
            NumberT distance,
            size_t bestTemplateIdx,
            const DynamicTimeWarping::MatchResult<NumberT>& matchResult,
            gsl::span<const DescriptorT> sequence,
            gsl::span<const double> timestamps,
            Session::Impl& sessionImpl)
        {
            Snapshot snapshot{};
            snapshot.timestamp = timestamp;
            snapshot.distance = distance;
            snapshot.imageStartIdx = matchResult.ImageStartIdx;
            snapshot.imageSize = matchResult.ImageSize;

            size_t startIdx = matchResult.ImageStartIdx;
            size_t endIdx = startIdx + matchResult.ImageSize;
            if (endIdx <= sequence.size())
            {
                snapshot.descriptors.assign(sequence.begin() + startIdx, sequence.begin() + endIdx);
            }

            if (startIdx < timestamps.size() && endIdx <= timestamps.size())
            {
                double startTime = timestamps[startIdx];
                double endTime = timestamps[endIdx - 1];
                if constexpr (std::is_same_v<DescriptorT, descriptor::Custom>)
                {
                    snapshot.inputSamples = sessionImpl.getSamplesInRange(startTime, endTime, ContractId<>::INVALID_ID);
                }
                else
                {
                    snapshot.inputSamples = sessionImpl.getSamplesInRange<DescriptorT>(startTime, endTime);
                }
            }

            m_pendingSnapshot = std::move(snapshot);
        }

        bool snapshotsCluster() const
        {
            if (m_snapshots.size() < 2)
            {
                return false;
            }

            auto absDist = [this](const DescriptorT& a, const DescriptorT& b) {
                if constexpr (descriptor::DescriptorTraits<DescriptorT>::HAS_ABSOLUTE_DISTANCE)
                {
                    return DescriptorT::AbsoluteDistance(a, b, m_tuning);
                }
                else
                {
                    return NumberT{ 0 };
                }
            };
            auto deltaDist = [this](const DescriptorT& a, const DescriptorT& a0, const DescriptorT& b, const DescriptorT& b0) {
                if constexpr (descriptor::DescriptorTraits<DescriptorT>::HAS_DELTA_DISTANCE)
                {
                    return DescriptorT::DeltaDistance(a, a0, b, b0, m_tuning);
                }
                else
                {
                    return NumberT{ 0 };
                }
            };

            NumberT totalDistance{};
            size_t pairCount{};

            for (size_t i = 0; i < m_snapshots.size(); ++i)
            {
                for (size_t j = i + 1; j < m_snapshots.size(); ++j)
                {
                    if (m_snapshots[i].descriptors.empty() || m_snapshots[j].descriptors.empty())
                    {
                        continue;
                    }

                    gsl::span<const DescriptorT> seqA{ m_snapshots[i].descriptors };
                    gsl::span<const DescriptorT> seqB{ m_snapshots[j].descriptors };
                    auto result = DynamicTimeWarping::Match(seqA, seqB, absDist, deltaDist);
                    if (result.Connections > 0)
                    {
                        totalDistance += result.MatchCost / result.Connections;
                        ++pairCount;
                    }
                }
            }

            if (pairCount == 0)
            {
                return false;
            }

            NumberT meanPairwiseDistance = totalDistance / pairCount;
            NumberT clusterThreshold = m_ambientEMA * m_settings.significanceMultiple;
            return meanPairwiseDistance < clusterThreshold;
        }
    };
}
