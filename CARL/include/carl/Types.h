/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <Eigen/Geometry>

namespace carl
{
    using NumberT = float;
    using VectorT = Eigen::Vector<NumberT, 3>;
    using QuaternionT = Eigen::Quaternion<NumberT>;
    using AngleAxisT = Eigen::AngleAxis<NumberT>;
    using TransformT = Eigen::Transform<NumberT, 3, Eigen::TransformTraits::Affine>;

    static inline const VectorT UNIT_SCALE{ 1, 1, 1 };

    namespace math
    {
        inline TransformT Lerp(const TransformT& a, const TransformT& b, NumberT t)
        {
            VectorT translation = a.translation() * (static_cast<NumberT>(1) - t) + b.translation() * t;
            QuaternionT rotation = QuaternionT{ a.rotation() }.slerp(t, QuaternionT{ b.rotation() });
            TransformT transform{};
            transform.fromPositionOrientationScale(translation, rotation, UNIT_SCALE);
            return transform;
        }

        inline TransformT LookTransform(const VectorT& forward, const VectorT& up, const VectorT& position)
        {
            const auto& yAxis = up;
            const auto& zAxis = forward;
            auto xAxis = yAxis.cross(zAxis);

            TransformT transform{};
            transform.matrix().col(0) = Eigen::Vector<NumberT, 4>(xAxis[0], xAxis[1], xAxis[2], static_cast<NumberT>(0));
            transform.matrix().col(1) = Eigen::Vector<NumberT, 4>(yAxis[0], yAxis[1], yAxis[2], static_cast<NumberT>(0));
            transform.matrix().col(2) = Eigen::Vector<NumberT, 4>(zAxis[0], zAxis[1], zAxis[2], static_cast<NumberT>(0));
            transform.matrix().col(3) = Eigen::Vector<NumberT, 4>(position[0], position[1], position[2], static_cast<NumberT>(1));
            return transform;
        }
    }
}
