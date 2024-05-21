/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "carl/Example.h"
#include "carl/Serialization.h"

namespace carl::action
{
    /// <summary>
    /// CARL's fundamental characterization of an action which can be recognized. An 
    /// action::Definition consists of a set of one or more Examples of the action
    /// being performed correctly, along with an ActionType indicating to what 
    /// category the define action belongs.
    /// </summary>
    class Definition
    {
    public:
        enum class ActionType : uint64_t
        {
            LeftHandPose = 0,
            LeftHandGesture = 1,
            RightHandPose = 2,
            RightHandGesture = 3,
            TwoHandGesture = 4,
            LeftControllerGesture = 5,
            RightControllerGesture = 6,
            TwoControllerGesture = 7,
        };

        Definition(ActionType actionType);
        Definition(Deserialization&);

        void serialize(Serialization&) const;

        void addExample(Example);
        void addCounterexample(Example);

        ActionType getDescriptorType() const
        {
            return m_actionType;
        }

        gsl::span<const Example> getExamples() const
        {
            return m_examples;
        }

        gsl::span<const Example> getCounterexamples() const
        {
            return m_counterexamples;
        }

        double DefaultSensitivity{ 1. };

    private:
        ActionType m_actionType{};
        std::vector<Example> m_examples{};
        std::vector<Example> m_counterexamples{};
    };
}
