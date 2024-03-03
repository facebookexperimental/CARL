/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/Session.h>

#include <arcana/threading/task.h>

#include "SessionImpl.h"

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
                resampling.emplace_back(InputSample::lerp(a, b, t));
            }
        }
    }

    Session::Impl::Impl(size_t samplesPerSecond, size_t maxActionDurationSeconds)
        : SessionImplBase{ samplesPerSecond * maxActionDurationSeconds }
        , frameDuration{ 1. / samplesPerSecond }
    {
    }

    Session::Impl& Session::Impl::getFromSession(Session& session)
    {
        return *session.m_impl;
    }

    void Session::Impl::addInputSample(const InputSample& inputSample)
    {
        // Send inputSample to the resampling logic. Only invoke the signal if the resampler updates the
        // input sequence.
        {
            std::scoped_lock lock{m_samplesMutex};
            appendSampleToResampling(inputSample, m_samples, frameDuration);
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

            if (m_processingSamples.size() > 1)
            {
                SignalHandlersT::apply_to_all([this](auto& callable) {
                    // The argument passed to the handlers here should be the NEW input samples -- only ones they
                    // haven't seen before PLUS the last one they HAVE seen before. This is a tricky and very
                    // tight coupling that exists between the session impl and the descriptor sequence providers,
                    // so it's okay as long as we contain it, but it should still and always be very clearly
                    // commented.
                    // TODO: Consider modifying the signal to make this l-value span unnecessary.
                    gsl::span<const InputSample> span{ m_processingSamples };
                    callable(span);
                });
            }
        });

        // TODO: Later, decouple this so that addInputSample is synchronous and resampling and processing
        // can happen on another thread.
    }

    Session::Session() : m_impl{ std::make_unique<Impl>() }
    {
    }

    Session::~Session()
    {
    }

    void Session::addInput(InputSample inputSample)
    {
        m_impl->addInputSample(inputSample);
    }

    /*arcana::background_dispatcher<128>& Session::processingScheduler()
    {
        return m_impl->processingScheduler();
    }*/

    arcana::manual_dispatcher<128>& Session::callbackScheduler()
    {
        return m_impl->callbackScheduler();
    }

    void Session::tickCallbacks(arcana::cancellation& token)
    {
        m_impl->tickCallbacks(token);
    }
}
