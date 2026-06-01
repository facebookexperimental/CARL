/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/ActionType.h>
#include <carl/Example.h>
#include <carl/Recording.h>

#include <gsl/span>

#include <vector>

namespace carl::action
{
    class DefinitionBuilder
    {
    public:
        /// Creates auto-trimmed Examples from multiple raw Recordings of
        /// the same unknown action. The recordings should each contain
        /// the action somewhere within them, potentially with leading and
        /// trailing excess. The algorithm discovers the common interesting
        /// sub-action across all recordings and returns trimmed Examples.
        ///
        /// @param actionType Determines which descriptor type to use.
        /// @param recordings Two or more raw recordings of the action.
        /// @param expectedActionDuration Optional hint (in seconds) for
        ///        the expected length of the action; 0 means unconstrained.
        /// @return A vector of trimmed Examples, one per input recording.
        static std::vector<Example> createExamplesFromRecordings(
            ActionType actionType,
            gsl::span<const Recording> recordings,
            double expectedActionDuration = 0.);
    };
}
