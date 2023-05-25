/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/Types.h>

namespace carl
{
    class EgocentricTemporalSpace
    {
    public:
        static TransformT getPose(VectorT origin, TransformT ego)
        {
            const auto& col0 = ego.matrix().col(0);
            const auto& col1 = ego.matrix().col(1);
            const VectorT up{ col1.x(), col1.y(), col1.z() };
            const VectorT right{ col0.x(), col0.y(), col0.z() };

            auto zAxis = (origin - ego.translation()).normalized();
            if (zAxis.dot(zAxis) < 0.00001f)
            {
                return TransformT::Identity(); // TODO DEBUG
            }

            auto yAccordingToEgoUp = zAxis.cross(up.cross(zAxis)).normalized();
            auto yAccordingToEgoRight = zAxis.cross(right).normalized();

            auto egoUpReliabilityFactor = 1.f - std::abs(zAxis.dot(up));
            auto egoRightReliabilityFactor = 1.f - std::abs(zAxis.dot(right));
            auto yAxis = (yAccordingToEgoUp * egoUpReliabilityFactor + yAccordingToEgoRight * egoRightReliabilityFactor).normalized();
            return math::LookTransform(zAxis, yAxis, origin);
        }
    };
}
