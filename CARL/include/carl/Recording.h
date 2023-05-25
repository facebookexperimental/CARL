/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "carl/InputSample.h"

#include <gsl/span>

namespace carl::action
{
    class Recording;

    /// <summary>
    /// A recording which is not yet completed. To create a new Recording,
    /// an InProgressRecording should first be created. Incoming InputSamples
    /// should then be added to the InProgressRecording, which can then be
    /// passed to (and subsumed into) a completed Recording.
    /// </summary>
    class InProgressRecording
    {
    public:
        InProgressRecording() = default;
        InProgressRecording(std::vector<InputSample> samples)
            : m_samples{ std::move(samples) }
        {
        }

        InProgressRecording(const InProgressRecording&) = delete;
        InProgressRecording(InProgressRecording&&) = default;
        InProgressRecording& operator=(const InProgressRecording&) = delete;
        InProgressRecording& operator=(InProgressRecording&&) = delete;

        void addSample(InputSample sample)
        {
            m_samples.emplace_back(std::move(sample));
        }

    private:
        friend class Recording;
        std::vector<InputSample> m_samples{};
    };

    /// <summary>
    /// A view on a Recording intended to make it easy and efficient to
    /// "scrub" through a Recording by timestamp. Note that RecordingInspector
    /// is strictly a view, with NO ownership of the underlying data, so a 
    /// RecordingInspector's lifespan must be strictly less than that of the
    /// data it was created to inspect.
    /// </summary>
    class RecordingInspector
    {
    public:
        RecordingInspector(gsl::span<const InputSample>);

        const InputSample& inspect(double timestamp);

        double startTimestamp() const
        {
            return m_samples.front().Timestamp;
        }

        double endTimestamp() const
        {
            return m_samples.back().Timestamp;
        }

    private:
        gsl::span<const InputSample> m_samples{};
        size_t m_idx{ 0 };
    };

    /// <summary>
    /// A time series of InputSamples. For a Recording to be valid, the 
    /// timestamps of all included InputSamples must be in reference to 
    /// the same point in time.
    /// </summary>
    class Recording
    {
    public:
        Recording(InProgressRecording);
        Recording(Deserialization&);

        Recording(const Recording&) = default;
        Recording(Recording&&) = default;
        Recording& operator=(const Recording&) = delete;
        Recording& operator=(Recording&&) = delete;

        void serialize(Serialization&) const;

        RecordingInspector getInspector() const
        {
            return{ m_samples };
        }

        gsl::span<const InputSample> getSamples() const
        {
            return m_samples;
        }

    private:
        std::vector<InputSample> m_samples{};
    };
}
