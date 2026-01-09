/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/ActionType.h>
#include <carl/Descriptor.h>

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
        static inline std::unique_ptr<action::Recognizer::Impl> createRecognizerForActionType(ActionType actionType, Session& session, const Definition& definition)
        {
            return internalCreateRecognizerForActionType<AssociationsT...>(actionType, session, definition);
        }

    private:
        template<typename AssociationT>
        static inline std::unique_ptr<action::Recognizer::Impl> createRecognizer(Session& session, const Definition& definition)
        {
            if constexpr (AssociationT::ExpandExamples)
            {
                auto examples = expandExamples(definition.getExamples());
                auto counterexamples = expandExamples(definition.getCounterexamples());
                return std::make_unique<RecognizerImpl<AssociationT::Descriptor>>(session, examples, counterexamples, definition.DefaultSensitivity);
            }
            else
            {
                return std::make_unique<RecognizerImpl<AssociationT::Descriptor>>(session, definition.getExamples(), definition.getCounterexamples(), definition.DefaultSensitivity);
            }
        }

        template<typename T>
        static inline std::unique_ptr<action::Recognizer::Impl> internalCreateRecognizerForActionType(ActionType actionType, Session& session, const Definition& definition)
        {
            if (actionType == T::Action)
            {
                return createRecognizer<T>(session, definition);
            }
            else
            {
                return {};
            }
        }

        template<typename T, typename T1, typename... Ts>
        static inline std::unique_ptr<action::Recognizer::Impl> internalCreateRecognizerForActionType(ActionType actionType, Session& session, const Definition& definition)
        {
            if (actionType == T::Action)
            {
                return createRecognizer<T>(session, definition);
            }
            else
            {
                return internalCreateRecognizerForActionType<T1, Ts...>(actionType, session, definition);
            }
        }
    };
}
