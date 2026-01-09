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
        template<typename RecognizerFactoriesT>
        static inline std::unique_ptr<action::Recognizer::Impl> createRecognizerForActionType(ActionType actionType, Session& session, const Definition& definition, const RecognizerFactoriesT& factories)
        {
            return internalCreateRecognizerForActionType<RecognizerFactoriesT, AssociationsT...>(actionType, session, definition, factories);
        }

    private:
        template<typename RecognizerFactoriesT, typename AssociationT>
        static inline std::unique_ptr<action::Recognizer::Impl> createRecognizer(Session& session, const Definition& definition, const RecognizerFactoriesT& factories)
        {
            auto& factory = factories.at(ContractId<typename AssociationT::Descriptor>::value());
            return factory(session, definition);
        }

        template<typename RecognizerFactoriesT, typename T>
        static inline std::unique_ptr<action::Recognizer::Impl> internalCreateRecognizerForActionType(ActionType actionType, Session& session, const Definition& definition, const RecognizerFactoriesT& factories)
        {
            if (actionType == T::Action)
            {
                return createRecognizer<RecognizerFactoriesT, T>(session, definition, factories);
            }
            else
            {
                return {};
            }
        }

        template<typename RecognizerFactoriesT, typename T, typename T1, typename... Ts>
        static inline std::unique_ptr<action::Recognizer::Impl> internalCreateRecognizerForActionType(ActionType actionType, Session& session, const Definition& definition, const RecognizerFactoriesT& factories)
        {
            if (actionType == T::Action)
            {
                return createRecognizer<RecognizerFactoriesT, T>(session, definition, factories);
            }
            else
            {
                return internalCreateRecognizerForActionType<RecognizerFactoriesT, T1, Ts...>(actionType, session, definition, factories);
            }
        }
    };
}
