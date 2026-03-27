/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

namespace carl::descriptor
{
    template<typename T>
    class DescriptorTraits
    {
        template<typename C>
        static char absoluteDistanceTest(decltype(&C::AbsoluteDistance));

        template<typename C>
        static char (&absoluteDistanceTest(...))[2];

        template<typename C>
        static char deltaDistanceTest(decltype(&C::DeltaDistance));

        template<typename C>
        static char (&deltaDistanceTest(...))[2];

    public:
        static constexpr bool HAS_ABSOLUTE_DISTANCE = sizeof(absoluteDistanceTest<T>(0)) == sizeof(char);
        static constexpr bool HAS_DELTA_DISTANCE = sizeof(deltaDistanceTest<T>(0)) == sizeof(char);
    };
}
