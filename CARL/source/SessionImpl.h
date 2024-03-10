/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "Descriptor.h"

#include <arcana/threading/dispatcher.h>

#include <carl/Example.h>
#include <carl/InputSample.h>
#include <carl/Recording.h>
#include <carl/Session.h>
#include <carl/Signaling.h>

#include <arcana/threading/task.h>

namespace carl
{
    template<typename DescriptorT, typename = std::enable_if_t<std::is_trivially_copyable_v<DescriptorT>>>
    class DescriptorSequence
    {
        using WeakTableT = arcana::weak_table<typename Signal<gsl::span<const DescriptorT>>::HandlerT>;

    public:
        using DescriptorSignalT = Signal<gsl::span<const DescriptorT>>;

        class Provider : private WeakTableT, public DescriptorSignalT
        {
        public:
            Provider(Signal<const InputSample&>& signal, size_t sequenceLength)
                : DescriptorSignalT{ *static_cast<WeakTableT*>(this) }
                , m_ticket{ signal.addHandler([this](auto samples) { handleInputSample(samples); }) }
                , m_sequenceLength{ sequenceLength }
            {
            }

            void handleInputSample(const InputSample& sample)
            {
                /*for (size_t idx = 1; idx < samples.size(); ++idx)
                {
                    auto descriptor = DescriptorT::TryCreate(samples[idx], samples[idx - 1]);
                    if (descriptor.has_value())
                    {
                        m_sequence.emplace_back(std::move(descriptor.value()));
                    }
                }*/
                size_t priorSequenceSize = m_sequence.size();
                descriptor::extendSequence(sample, m_sequence, m_mostRecentSample, DescriptorT::DEFAULT_TUNING);
                bool descriptorsAdded = m_sequence.size() > priorSequenceSize;

                if (m_sequence.size() > m_sequenceLength)
                {
                    // Note: This does dangerous direct memory manipulation to avoid calling a bunch of
                    // range calls. This is allowable because of the requirement that descriptors be
                    // trivially copyable, but the behavior has the potential to be subtle and should
                    // be treated with caution.
                    m_buffer.resize(m_sequenceLength);
                    std::memcpy(m_buffer.data(), m_sequence.data() + (m_sequence.size() - m_sequenceLength), sizeof(DescriptorT) * m_sequenceLength);
                    m_buffer.swap(m_sequence);
                }

                if (descriptorsAdded)
                {
                    auto& handlers = *static_cast<WeakTableT*>(this);
                    handlers.apply_to_all([this](auto& callable) {
                        gsl::span<const DescriptorT> span{ m_sequence };
                        callable(span);
                    });
                }
            }

        private:
            const Signal<const InputSample&>::TicketT m_ticket;
            InputSample m_mostRecentSample{};
            std::vector<DescriptorT> m_sequence{};
            std::vector<DescriptorT> m_buffer{};
            const size_t m_sequenceLength{};
        };
    };

    template <typename... DescriptorTs>
    class SessionImplBase
        : public arcana::weak_table<typename Signal<const InputSample&>::HandlerT>
        , public Signal<const InputSample&>
        , protected DescriptorSequence<DescriptorTs>::Provider...
    {
    public:
        using SignalHandlersT = arcana::weak_table<typename Signal<const InputSample&>::HandlerT>;

        SessionImplBase(size_t sequenceLength)
            : Signal<const InputSample&>{ *static_cast<SignalHandlersT*>(this) }
            , DescriptorSequence<DescriptorTs>::Provider{ *static_cast<Signal<const InputSample&>*>(this), sequenceLength }...
        {
        }
    };

    class Session::Impl
        : public SessionImplBase<
        descriptor::HandPose<descriptor::Handedness::LeftHanded>,
        descriptor::HandGesture<descriptor::Handedness::LeftHanded>,
        descriptor::HandPose<descriptor::Handedness::RightHanded>,
        descriptor::HandGesture<descriptor::Handedness::RightHanded>,
        descriptor::TwoHandGesture>
    {
    public:
        Impl(size_t samplesPerSecond, size_t maxActionDurationSeconds, bool singleThreaded);

        static Session::Impl& getFromSession(Session& session);

        void addInputSample(const InputSample& inputSample);

        auto& processingScheduler()
        {
            return m_processingScheduler;
        }

        auto& callbackScheduler()
        {
            return m_callbackScheduler;
        }

        void tickCallbacks(arcana::cancellation& token)
        {
            m_callbackDispatcher.tick(token);
        }

        const double frameDuration{};

        template <typename DescriptorT>
        auto addHandler(std::function<void(gsl::span<const DescriptorT>)> handler)
        {
            return typename DescriptorSequence<DescriptorT>::Provider::addHandler(std::move(handler));
        }

        void setLogger(std::function<void(std::string)> logger)
        {
            std::scoped_lock lock{ m_loggerMutex };
            m_logger = std::move(logger);
        }

    private:
        arcana::manual_dispatcher<256> m_callbackDispatcher{};
        arcana::background_dispatcher<256> m_processingDispatcher{};
        SchedulerT m_callbackScheduler{};
        SchedulerT m_processingScheduler{};
        std::vector<InputSample> m_samples{};
        std::vector<InputSample> m_processingSamples{};
        std::mutex m_samplesMutex{};
        std::function<void(std::string)> m_logger{};
        std::mutex m_loggerMutex{};
    };
}
