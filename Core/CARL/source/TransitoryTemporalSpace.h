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
    class TransitoryTemporalSpace
    {
    public:
        static TransformT getPose(VectorT origin, VectorT priorOrigin, TransformT ego)
        {
            auto zAxis = (origin - priorOrigin).normalized();
            auto yAxis = zAxis.cross(ego.translation() - origin).normalized();
            return math::LookTransform(zAxis, yAxis, origin);
        }
    };
}
