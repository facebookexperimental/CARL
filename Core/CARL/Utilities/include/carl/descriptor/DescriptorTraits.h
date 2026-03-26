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
        static char test(decltype(&C::DeltaDistance));

        template<typename C>
        static char (&test(...))[2];

    public:
        static constexpr bool IS_DELTA = sizeof(test<T>(0)) == sizeof(char);
        static constexpr bool HAS_DELTA = IS_DELTA;
    };
}
