/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <atomic>

namespace carl
{
    template<typename...>
    class ContractId;

    namespace internal
    {
        class ContractIdValue final
        {
        private:
            using IdT = size_t;

            template <typename...>
            friend class carl::ContractId;

            static size_t next()
            {
                static std::atomic<size_t> nextId{ 1 };
                return nextId.fetch_add(1);
            }
        };
    }

    template<>
    class ContractId<>
    {
    public:
        using IdT = internal::ContractIdValue::IdT;
        static inline constexpr IdT INVALID_ID{ 0 };
        static IdT reserveDynamicContractId()
        {
            return internal::ContractIdValue::next();
        }
    };

    template<typename T>
    class ContractId<T>
    {
    public:
        using IdT = internal::ContractIdValue::IdT;
        static inline constexpr IdT INVALID_ID{ 0 };
        static IdT value()
        {
            static auto value{ internal::ContractIdValue::next() };
            return value;
        }
    };
}