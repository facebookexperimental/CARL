/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/InputSample.h>

#include <vector>

namespace carl::action
{
    struct ActionAttemptEstimationSettings
    {
        double ambientAlpha{ 0.02 };
        double significanceMultiple{ 0.6 };
        double clusterWindowSeconds{ 5.0 };
        size_t clusterMinCount{ 3 };
        double snapshotExpirySeconds{ 10.0 };
        std::optional<double> precomputedSignificanceThreshold{};
    };

    struct AttemptCluster
    {
        std::vector<std::vector<InputSample>> snapshots{};
        double firstTimestamp{};
        double lastTimestamp{};
        size_t count{};
    };
}
