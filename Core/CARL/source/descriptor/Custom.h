/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "DescriptorUtils.h"

namespace carl::descriptor
{
    class Custom
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Custom";

    public:
        static constexpr std::array<NumberT, 32> DEFAULT_TUNING{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

        Custom(std::shared_ptr<void> descriptor, const internal::CustomActionTypeOperations& operations)
            : m_descriptor{ std::move(descriptor) }
            , m_operations{ operations }
        {
        }

        Custom(const Custom& other)
            : m_descriptor{ other.m_descriptor }
            , m_operations{ other.m_operations }
        {
        }

        Custom& operator=(const Custom& other)
        {
            m_descriptor = other.m_descriptor;
            assert(&m_operations == &other.m_operations);
            return *this;
        }

        static NumberT Distance(const Custom& a, const Custom& a0, const Custom& b, const Custom& b0, gsl::span<const NumberT> tuning)
        {
            return a.m_operations.Distance(a.m_descriptor, a0.m_descriptor, b.m_descriptor, b0.m_descriptor, tuning);
        }

        static Custom Lerp(const Custom& a, const Custom& b, NumberT t)
        {
            assert(&a.m_operations == &b.m_operations);
            auto desc = a.m_operations.Lerp(a.m_descriptor, b.m_descriptor, t);
            return{ std::move(desc), a.m_operations };
        }

    private:
        std::shared_ptr<void> m_descriptor{};
        const internal::CustomActionTypeOperations& m_operations;
    };
}
