/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "carl/InputSample.h"

#include <memory>

namespace carl
{
    /// <summary>
    /// CARL's characterization of a timespan and usage in which actions can be recognized.
    /// </summary>
    class Session {
    public:
        class Impl;

        Session();
        ~Session();

        void addInput(InputSample);

    private:
        friend class Impl;
        std::unique_ptr<Impl> m_impl{};
    };
}
