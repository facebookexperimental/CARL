/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "carl/InputSample.h"

#include <arcana/threading/dispatcher.h>

#include <memory>

namespace carl
{
    /// <summary>
    /// CARL's characterization of a timespan and usage in which actions can be recognized.
    /// </summary>
    class Session {
    public:
        class Impl;
        using SchedulerT = stdext::inplace_function<void(stdext::inplace_function<void(), 128>&&), 128>;

        Session(bool singleThreaded = false);
        ~Session();

        void addInput(InputSample);

        SchedulerT& callbackScheduler();
        SchedulerT& processingScheduler();
        void setLogger(std::function<void(std::string)> logger);
        void log(std::string message);
        void tickCallbacks(arcana::cancellation& token);

    private:
        friend class Impl;
        std::unique_ptr<Impl> m_impl{};
    };
}
