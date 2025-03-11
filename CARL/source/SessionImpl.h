/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "Descriptor.h"

#include <carl/Example.h>
#include <carl/InputSample.h>
#include <carl/Recording.h>
#include <carl/Session.h>
#include <carl/Signaling.h>

#include <arcana/threading/dispatcher.h>
#include <arcana/threading/task.h>

namespace carl
{
    template<typename DescriptorT, typename = std::enable_if_t<std::is_trivially_copyable_v<DescriptorT>>>
    class DescriptorSequence
    {
        using WeakTableT = arcana::weak_table<typename Signal<gsl::span<const DescriptorT>, size_t>::HandlerT>;

    public:
        using DescriptorSignalT = Signal<gsl::span<const DescriptorT>, size_t>;

        class Provider : private WeakTableT, public DescriptorSignalT
        {
        public:
            Provider(Signal<gsl::span<const InputSample>>& signal)
                : DescriptorSignalT{ *static_cast<WeakTableT*>(this) }
                , m_ticket{ signal.addHandler([this](auto samples) { handleInputSamples(samples); }) }
                , m_sequenceLength{ 10 }
            {
            }

            void handleInputSamples(gsl::span<const InputSample> samples)
            {
                size_t descriptorsAddedCount = 0;
                for (const auto& sample : samples)
                {
                    descriptorsAddedCount += handleInputSample(sample);
                }

                if (descriptorsAddedCount > 0)
                {
                    auto& handlers = *static_cast<WeakTableT*>(this);
                    gsl::span<const DescriptorT> span{ m_sequence };
                    handlers.apply_to_all([this, span, descriptorsAddedCount](auto& callable) mutable {
                        callable(span, descriptorsAddedCount);
                    });
                }
            }

            size_t handleInputSample(const InputSample& sample)
            {
                size_t priorSequenceSize = m_sequence.size();
                descriptor::extendSequence(sample, m_sequence, m_mostRecentSample, DescriptorT::DEFAULT_TUNING);
                size_t descriptorsAddedCount = m_sequence.size() - priorSequenceSize;

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

                return descriptorsAddedCount;
            }

            void supportSequenceOfLength(size_t length)
            {
                m_sequenceLength = std::max(m_sequenceLength, length);
            }

        private:
            const Signal<gsl::span<const InputSample>>::TicketT m_ticket;
            InputSample m_mostRecentSample{};
            std::vector<DescriptorT> m_sequence{};
            std::vector<DescriptorT> m_buffer{};
            size_t m_sequenceLength{};
        };
    };

    template <typename... DescriptorTs>
    class SessionImplBase
        : public arcana::weak_table<typename Signal<gsl::span<const InputSample>>::HandlerT>
        , public Signal<gsl::span<const InputSample>>
        , protected DescriptorSequence<DescriptorTs>::Provider...
    {
    public:
        using SignalHandlersT = arcana::weak_table<typename Signal<gsl::span<const InputSample>>::HandlerT>;

        SessionImplBase()
            : Signal<gsl::span<const InputSample>>{ *static_cast<SignalHandlersT*>(this) }
            , DescriptorSequence<DescriptorTs>::Provider{ *static_cast<Signal<gsl::span<const InputSample>>*>(this) }...
        {
        }
    };

    class Session::Impl
        : public SessionImplBase<
        descriptor::HandPose<descriptor::Handedness::LeftHanded>,
        descriptor::HandGesture<descriptor::Handedness::LeftHanded>,
        descriptor::HandPose<descriptor::Handedness::RightHanded>,
        descriptor::HandGesture<descriptor::Handedness::RightHanded>,
        descriptor::TwoHandGesture,
        descriptor::ControllerGesture<descriptor::Handedness::LeftHanded>,
        descriptor::ControllerGesture<descriptor::Handedness::RightHanded>,
        descriptor::TwoControllerGesture,
        descriptor::WristTrajectory<descriptor::Handedness::LeftHanded>,
        descriptor::WristTrajectory<descriptor::Handedness::RightHanded>,
        descriptor::HandShape<descriptor::Handedness::LeftHanded>,
        descriptor::HandShape<descriptor::Handedness::RightHanded>>
    {
    public:
        Impl(bool singleThreaded);

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

        template <typename DescriptorT>
        auto addHandler(std::function<void(gsl::span<const DescriptorT>, size_t)> handler)
        {
            return DescriptorSequence<DescriptorT>::Provider::addHandler(std::move(handler));
        }

        template<typename DescriptorT>
        void supportSequenceOfLength(size_t length)
        {
            DescriptorSequence<DescriptorT>::Provider::supportSequenceOfLength(length);
        }

        void setLogger(std::function<void(std::string)> logger)
        {
            std::scoped_lock lock{ m_loggerMutex };
            m_logger = std::move(logger);
        }

        void log(std::string message)
        {
            std::scoped_lock lock{ m_loggerMutex };
            m_logger(std::move(message));
        }

    private:
        arcana::manual_dispatcher<256> m_callbackDispatcher{};
        std::optional<arcana::background_dispatcher<256>> m_processingDispatcher{};
        SchedulerT m_callbackScheduler{};
        SchedulerT m_processingScheduler{};
        std::vector<InputSample> m_samples{};
        std::vector<InputSample> m_processingSamples{};
        std::mutex m_samplesMutex{};
        std::function<void(std::string)> m_logger{};
        std::mutex m_loggerMutex{};
    };
}
