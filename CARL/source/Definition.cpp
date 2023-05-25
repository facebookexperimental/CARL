/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/Definition.h>

namespace carl::action
{
    Definition::Definition(ActionType descriptorType)
        : m_actionType{ descriptorType }
    {
    }

    Definition::Definition(Deserialization& deserialization)
    {
        deserialization >> m_actionType;
        uint64_t count;
        deserialization >> count;
        m_examples.reserve(count);
        for (uint64_t idx = 0; idx < count; ++idx)
        {
            m_examples.emplace_back(deserialization);
        }
        deserialization >> DefaultSensitivity;
    }

    void Definition::serialize(Serialization& serialization) const
    {
        serialization << m_actionType;
        serialization << static_cast<uint64_t>(m_examples.size());
        for (const auto& example : m_examples)
        {
            example.serialize(serialization);
        }
        serialization << DefaultSensitivity;
    }

    void Definition::addExample(Example example)
    {
        m_examples.emplace_back(std::move(example));
    }
}
