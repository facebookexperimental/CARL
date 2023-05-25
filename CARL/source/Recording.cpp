/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/Recording.h>

namespace carl::action {
    RecordingInspector::RecordingInspector(gsl::span<const InputSample> samples)
        : m_samples{ samples }
    {
    }

    const InputSample& RecordingInspector::inspect(double timestamp)
    {
        // Extremely naively walk the samples vector.
        while (true)
        {
            if (m_samples[m_idx].Timestamp <= timestamp)
            {
                // If incrementing would go out of bounds, halt.
                if (m_idx + 1 >= m_samples.size())
                {
                    break;
                }

                // If incrementing would move past the target timestamp, halt.
                if (m_samples[m_idx + 1].Timestamp > timestamp)
                {
                    break;
                }

                ++m_idx;
            }
            else
            {
                // If decrementing would go out of bounds, halt.
                if (m_idx == 0)
                {
                    break;
                }

                --m_idx;
            }
        }

        return m_samples[m_idx];
    }

    Recording::Recording(InProgressRecording recording) 
        : m_samples{ std::move(recording.m_samples) }
    {
        std::sort(m_samples.begin(), m_samples.end(), [](const auto& a, const auto& b) {
            return a.Timestamp < b.Timestamp;
        });
    }

    Recording::Recording(Deserialization& deserialization)
    {
        uint64_t count;
        deserialization >> count;
        m_samples.reserve(count);
        for (uint64_t idx = 0; idx < count; ++idx)
        {
            m_samples.emplace_back(deserialization);
        }
    }

    void Recording::serialize(Serialization& serialization) const
    {
        serialization << static_cast<uint64_t>(m_samples.size());
        for (const auto& sample : m_samples)
        {
            sample.serialize(serialization);
        }
    }
}
