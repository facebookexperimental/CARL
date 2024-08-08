/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/Session.h>

#include "SessionImpl.h"
#include "Descriptor.h"

#include <arcana/threading/task.h>

namespace carl
{
    namespace
    {
        // Naive incremental linear resampling.
        void appendSampleToResampling(const InputSample& sample, std::vector<InputSample>& resampling, double frameDuration)
        {
            if (resampling.empty())
            {
                resampling.emplace_back(sample);
            }

            while (resampling.back().Timestamp + frameDuration <= sample.Timestamp)
            {
                const auto& a = resampling.back();
                const auto& b = sample;
                double timestamp = a.Timestamp + frameDuration;
                double t = (timestamp - a.Timestamp) / (b.Timestamp - a.Timestamp);
                resampling.emplace_back(InputSample::Lerp(a, b, t));
            }
        }
    }

    Session::Impl::Impl(bool singleThreaded)
        : m_callbackScheduler{ [this](auto&& work) { m_callbackDispatcher(std::forward<std::remove_reference_t<decltype(work)>>(work)); }}
        , m_processingScheduler{ [this, singleThreaded]() -> SchedulerT {
            if (singleThreaded)
            {
                return [](auto&& work) { 
                    arcana::inline_scheduler(std::forward<std::remove_reference_t<decltype(work)>>(work));
                };
            }
            else
            {
                m_processingDispatcher.emplace();
                return [&dispatcher = m_processingDispatcher.value()](auto&& work) {
                    dispatcher(std::forward<std::remove_reference_t<decltype(work)>>(work));
                };
            }
        }()}
    {
    }

    Session::Impl& Session::Impl::getFromSession(Session& session)
    {
        return *session.m_impl;
    }

    void Session::Impl::addInputSample(const InputSample& inputSample)
    {
        {
            std::scoped_lock lock{m_samplesMutex};
            m_samples.push_back(inputSample);
        }
        arcana::make_task(processingScheduler(), arcana::cancellation::none(), [this]() {
            {
                std::scoped_lock lock{m_samplesMutex};
                if (m_samples.size() > 1) {
                    m_processingSamples.swap(m_samples);
                    m_samples.resize(1);
                    m_samples[0] = m_processingSamples.back();
                }
            }

            gsl::span<const InputSample> span{ m_processingSamples };
            SignalHandlersT::apply_to_all([span](auto& callable) mutable {
                callable(span);
            });
        }).then(arcana::inline_scheduler, arcana::cancellation::none(), [this](arcana::expected<void, std::exception_ptr> expected) {
            if (expected.has_error())
            {
                try
                {
                    std::rethrow_exception(expected.error());
                }
                catch (std::exception& e)
                {
                    arcana::make_task(m_callbackScheduler, arcana::cancellation::none(), [this, message = std::string{ e.what() }]() {
                        std::scoped_lock lock{ m_loggerMutex };
                        m_logger(message.c_str());
                    });
                }
            }
        });
    }

    Session::Session(bool singleThreaded) 
        : m_impl{ std::make_unique<Impl>(singleThreaded) }
    {
    }

    Session::~Session()
    {
    }

    void Session::addInput(const InputSample& inputSample)
    {
        m_impl->addInputSample(inputSample);
    }

    Session::SchedulerT& Session::processingScheduler()
    {
        return m_impl->processingScheduler();
    }

    Session::SchedulerT& Session::callbackScheduler()
    {
        return m_impl->callbackScheduler();
    }

    void Session::setLogger(std::function<void(std::string)> logger)
    {
        m_impl->setLogger(std::move(logger));
    }

    void Session::log(std::string message)
    {
        m_impl->log(std::move(message));
    }

    void Session::tickCallbacks(arcana::cancellation& token)
    {
        m_impl->tickCallbacks(token);
    }
}
