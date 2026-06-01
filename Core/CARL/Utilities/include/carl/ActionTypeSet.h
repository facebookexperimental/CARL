/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/ActionType.h>
#include <carl/ContractId.h>
#include <carl/descriptor/SequenceOperations.h>

namespace carl::action
{
    template <ActionType ActionT, typename DescriptorT, bool expandExamples = false>
    struct ActionTypeDefinition
    {
        static inline constexpr ActionType Action{ ActionT };
        using Descriptor = DescriptorT;
        static inline constexpr bool ExpandExamples{ expandExamples };
    };

    template<typename... AssociationsT>
    struct ActionTypeSet
    {
        template<typename CallableT>
        static inline auto invokeByActionType(ActionType actionType, CallableT&& callable)
        {
            return internalInvokeByActionType<AssociationsT...>(actionType, std::forward<CallableT>(callable));
        }

    private:
        template<typename T, typename CallableT>
        static inline auto internalInvokeByActionType(ActionType actionType, CallableT&& callable)
            -> decltype(callable(static_cast<T*>(nullptr)))
        {
            if (actionType == T::Action)
            {
                return callable(static_cast<T*>(nullptr));
            }
            return {};
        }

        template<typename T, typename T1, typename... Ts, typename CallableT>
        static inline auto internalInvokeByActionType(ActionType actionType, CallableT&& callable)
            -> decltype(callable(static_cast<T*>(nullptr)))
        {
            if (actionType == T::Action)
            {
                return callable(static_cast<T*>(nullptr));
            }
            return internalInvokeByActionType<T1, Ts...>(actionType, std::forward<CallableT>(callable));
        }
    };
}
