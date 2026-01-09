/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/ContractId.h>
#include <arcana/functional/inplace_function.h>

#include <unordered_map>

namespace carl
{
    class TypedCollection
    {
    public:
        template<typename T, typename... Ts>
        T& getOrCreateInstanceForContractId(ContractId<>::IdT id, Ts&&... args)
        {
            auto found = m_contractIdToInstance.find(id);
            if (found == m_contractIdToInstance.end())
            {
                auto outerPtr = std::make_unique<T>(std::forward<Ts>(args)...);
                auto inserted = m_contractIdToInstance.try_emplace(id, [ptr{ std::move(outerPtr) }]() -> void* {
                    return ptr.get();
                    });
                return *reinterpret_cast<T*>(inserted.first->second());
            }
            else
            {
                return *reinterpret_cast<T*>(found->second());
            }
        }

        template<typename T, typename... Ts>
        T& getOrCreateInstance(Ts&&... args)
        {
            return getOrCreateInstanceForContractId<T>(ContractId<T>::value(), std::forward<Ts>(args)...);
        }

    private:
        std::unordered_map<size_t, stdext::inplace_function<void* (), stdext::InplaceFunctionDefaultCapacity, alignof(std::max_align_t), false>> m_contractIdToInstance{};
    };
}
