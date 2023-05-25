/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "carl/Recording.h"
#include "carl/Serialization.h"

namespace carl::action
{
    /// <summary>
    /// An Example is essentially a Recording in which it is known that something
    /// interesting (i.e., an action) transpired in a timespan bounded by 
    /// startTimestamp and endTimestamp. In some sense, an Example can be thought
    /// of as a "trimmed" recording.
    /// </summary>
    class Example
    {
    public:
        Example(Recording recording, double startTimestamp, double endTimestamp);
        Example(Deserialization&);

        void serialize(Serialization&) const;

        const Recording& getRecording() const
        {
            return m_recording;
        }

        double getStartTimestamp() const
        {
            return m_startTimestamp;
        }

        double getEndTimestamp() const
        {
            return m_endTimestamp;
        }

    private:
        Recording m_recording;
        double m_startTimestamp{};
        double m_endTimestamp{};
    };
}
