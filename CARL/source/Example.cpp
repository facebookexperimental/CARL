/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/Example.h>

namespace carl::action
{
    Example::Example(Recording recording, double startTimestamp, double endTimestamp)
        : m_recording{ std::move(recording) }
        , m_startTimestamp{ startTimestamp }
        , m_endTimestamp{ endTimestamp }
    {
    }

    Example::Example(Deserialization& deserialization)
        : m_recording{ deserialization }
    {
        deserialization >> m_startTimestamp;
        deserialization >> m_endTimestamp;
    }

    void Example::serialize(Serialization& serialization) const
    {
        m_recording.serialize(serialization);
        serialization << m_startTimestamp;
        serialization << m_endTimestamp;
    }
}
