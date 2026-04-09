/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/descriptor/Custom.h>
#include <carl/descriptor/SequenceOperations.h>

#include <carl/Example.h>
#include <carl/InputSample.h>
#include <carl/Recording.h>
#include <carl/RingBuffer.h>
#include <carl/Session.h>
#include <carl/Signaling.h>
#include <carl/TypedCollection.h>

#include <arcana/threading/dispatcher.h>
#include <arcana/threading/task.h>

namespace carl
{
    template<typename DescriptorT>
    class DescriptorSequence
    {
        using WeakTableT = arcana::weak_table<typename Signal<gsl::span<const DescriptorT>, size_t>::HandlerT>;

    public:
        using DescriptorSignalT = Signal<gsl::span<const DescriptorT>, size_t>;

        class Provider : private WeakTableT, public DescriptorSignalT
        {
        public:
            Provider(Signal<gsl::span<const InputSample>>& signal, internal::CustomActionTypeOperations* operations = nullptr)
                : DescriptorSignalT{ *static_cast<WeakTableT*>(this) }
            , m_ticket{ signal.addHandler([this](auto samples) { handleInputSamples(samples); }) }
                , m_sequenceLength{ 10 }
                , m_operations{ operations }
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
                if (m_retainedSamples.has_value())
                {
                    m_retainedSamples->push_back(sample);
                }

                size_t priorSequenceSize = m_sequence.size();
                if constexpr (std::is_same_v<DescriptorT, descriptor::Custom>)
                {
                    auto tryCreate = [this](const InputSample& current, const InputSample& prior) -> std::optional<descriptor::Custom> {
                        auto desc = m_operations->TryCreate(current, prior);
                        if (desc.has_value())
                        {
                            return descriptor::Custom{ std::move(desc.value()), *m_operations };
                        }
                        else
                        {
                            return{};
                        }
                        };
                    descriptor::extendSequence(sample, m_sequence, m_mostRecentSample, DescriptorT::DEFAULT_TUNING, tryCreate, m_timestamps);
                }
                else
                {
                    constexpr auto tryCreate = [](const InputSample& newSample, const InputSample& priorSample) { return DescriptorT::TryCreate(newSample, priorSample); };
                    descriptor::extendSequence(sample, m_sequence, m_mostRecentSample, DescriptorT::DEFAULT_TUNING, tryCreate, m_timestamps);
                }
                size_t descriptorsAddedCount = m_sequence.size() - priorSequenceSize;

                if (m_sequence.size() > m_sequenceLength)
                {
                    size_t excess = m_sequence.size() - m_sequenceLength;
                    if constexpr (std::is_trivially_copyable_v<DescriptorT>)
                    {
                        // Note: This does dangerous direct memory manipulation to avoid calling a bunch of
                        // range calls. This is allowable because of the requirement that descriptors be
                        // trivially copyable, but the behavior has the potential to be subtle and should
                        // be treated with caution.
                        m_buffer.resize(m_sequenceLength);
                        std::memcpy(m_buffer.data(), m_sequence.data() + excess, sizeof(DescriptorT) * m_sequenceLength);
                        m_buffer.swap(m_sequence);
                    }
                    else
                    {
                        m_sequence.erase(m_sequence.begin(), m_sequence.begin() + excess);
                    }

                    m_timestampBuffer.resize(m_sequenceLength);
                    std::memcpy(m_timestampBuffer.data(), m_timestamps.data() + excess, sizeof(double) * m_sequenceLength);
                    m_timestampBuffer.swap(m_timestamps);
                }

                return descriptorsAddedCount;
            }

            void supportSequenceOfLength(size_t length)
            {
                m_sequenceLength = std::max(m_sequenceLength, length);
            }

            gsl::span<const double> getTimestamps() const
            {
                return m_timestamps;
            }

            void enableSampleRetention(size_t capacity)
            {
                if (!m_retainedSamples.has_value())
                {
                    m_retainedSamples.emplace(capacity);
                }
            }

            const std::optional<RingBuffer<InputSample>>& getRetainedSamples() const
            {
                return m_retainedSamples;
            }

        private:
            const Signal<gsl::span<const InputSample>>::TicketT m_ticket;
            InputSample m_mostRecentSample{};
            std::vector<DescriptorT> m_sequence{};
            std::vector<DescriptorT> m_buffer{};
            std::vector<double> m_timestamps{};
            std::vector<double> m_timestampBuffer{};
            size_t m_sequenceLength{};
            const internal::CustomActionTypeOperations* m_operations{};
            std::optional<RingBuffer<InputSample>> m_retainedSamples{};
        };
    };

    class SessionImplBase
        : public arcana::weak_table<typename Signal<gsl::span<const InputSample>>::HandlerT>
        , public Signal<gsl::span<const InputSample>>
    {
    public:
        using SignalHandlersT = arcana::weak_table<typename Signal<gsl::span<const InputSample>>::HandlerT>;

        SessionImplBase()
            : Signal<gsl::span<const InputSample>>{ *static_cast<SignalHandlersT*>(this) }
        {
        }
    };

    class Session::Impl : public SessionImplBase
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

        ContractId<>::IdT enableCustomActionType(internal::CustomActionTypeOperations ops);

        const internal::CustomActionTypeOperations& getCustomActionTypeOperations(ContractId<>::IdT id) const
        {
            return *m_customActionIdsToOperations.at(id);
        }

        template<typename DescriptorT, typename = std::enable_if_t<!std::is_same_v<DescriptorT, descriptor::Custom>>>
        auto addHandler(std::function<void(gsl::span<const DescriptorT>, size_t)> handler)
        {
            carl::Signal<gsl::span<const InputSample>>& signal{ *this };
            auto& provider = m_dynamicSequenceProviders.getOrCreateInstance<typename DescriptorSequence<DescriptorT>::Provider>(signal);
            if (m_sampleRetentionCapacity.has_value())
            {
                provider.enableSampleRetention(*m_sampleRetentionCapacity);
            }
            return provider.addHandler(std::move(handler));
        }

        auto addHandler(std::function<void(gsl::span<const descriptor::Custom>, size_t)> handler, ContractId<>::IdT contractId)
        {
            carl::Signal<gsl::span<const InputSample>>& signal{ *this };
            auto* operations = m_customActionIdsToOperations.at(contractId).get();
            auto& provider = m_dynamicSequenceProviders.getOrCreateInstanceForContractId<typename DescriptorSequence<descriptor::Custom>::Provider>(contractId, signal, operations);
            if (m_sampleRetentionCapacity.has_value())
            {
                provider.enableSampleRetention(*m_sampleRetentionCapacity);
            }
            return provider.addHandler(std::move(handler));
        }

        template<typename DescriptorT, typename = std::enable_if_t<!std::is_same_v<DescriptorT, descriptor::Custom>>>
        void supportSequenceOfLength(size_t length)
        {
            carl::Signal<gsl::span<const InputSample>>& signal{ *this };
            auto& provider = m_dynamicSequenceProviders.getOrCreateInstance<typename DescriptorSequence<DescriptorT>::Provider>(signal);
            return provider.supportSequenceOfLength(length);
        }

        void supportSequenceOfLength(size_t length, ContractId<>::IdT contractId)
        {
            carl::Signal<gsl::span<const InputSample>>& signal{ *this };
            auto* operations = m_customActionIdsToOperations.at(contractId).get();
            auto& provider = m_dynamicSequenceProviders.getOrCreateInstanceForContractId<typename DescriptorSequence<descriptor::Custom>::Provider>(contractId, signal, operations);
            return provider.supportSequenceOfLength(length);
        }

        template<typename DescriptorT, typename = std::enable_if_t<!std::is_same_v<DescriptorT, descriptor::Custom>>>
        gsl::span<const double> getProviderTimestamps()
        {
            carl::Signal<gsl::span<const InputSample>>& signal{ *this };
            auto& provider = m_dynamicSequenceProviders.getOrCreateInstance<typename DescriptorSequence<DescriptorT>::Provider>(signal);
            return provider.getTimestamps();
        }

        gsl::span<const double> getProviderTimestamps(ContractId<>::IdT contractId)
        {
            carl::Signal<gsl::span<const InputSample>>& signal{ *this };
            auto* operations = m_customActionIdsToOperations.at(contractId).get();
            auto& provider = m_dynamicSequenceProviders.getOrCreateInstanceForContractId<typename DescriptorSequence<descriptor::Custom>::Provider>(contractId, signal, operations);
            return provider.getTimestamps();
        }

        template<typename DescriptorT, typename = std::enable_if_t<!std::is_same_v<DescriptorT, descriptor::Custom>>>
        void enableSampleRetention(size_t capacity)
        {
            carl::Signal<gsl::span<const InputSample>>& signal{ *this };
            auto& provider = m_dynamicSequenceProviders.getOrCreateInstance<typename DescriptorSequence<DescriptorT>::Provider>(signal);
            provider.enableSampleRetention(capacity);
        }

        void enableSampleRetention(size_t capacity, ContractId<>::IdT contractId)
        {
            carl::Signal<gsl::span<const InputSample>>& signal{ *this };
            auto* operations = m_customActionIdsToOperations.at(contractId).get();
            auto& provider = m_dynamicSequenceProviders.getOrCreateInstanceForContractId<typename DescriptorSequence<descriptor::Custom>::Provider>(contractId, signal, operations);
            provider.enableSampleRetention(capacity);
        }

        template<typename DescriptorT, typename = std::enable_if_t<!std::is_same_v<DescriptorT, descriptor::Custom>>>
        std::vector<InputSample> getSamplesInRange(double startTime, double endTime)
        {
            carl::Signal<gsl::span<const InputSample>>& signal{ *this };
            auto& provider = m_dynamicSequenceProviders.getOrCreateInstance<typename DescriptorSequence<DescriptorT>::Provider>(signal);
            return extractSamplesFromRetained(provider.getRetainedSamples(), startTime, endTime);
        }

        std::vector<InputSample> getSamplesInRange(double startTime, double endTime, ContractId<>::IdT contractId)
        {
            carl::Signal<gsl::span<const InputSample>>& signal{ *this };
            auto* operations = m_customActionIdsToOperations.at(contractId).get();
            auto& provider = m_dynamicSequenceProviders.getOrCreateInstanceForContractId<typename DescriptorSequence<descriptor::Custom>::Provider>(contractId, signal, operations);
            return extractSamplesFromRetained(provider.getRetainedSamples(), startTime, endTime);
        }

        void setLogger(std::function<void(std::string)> logger)
        {
            std::scoped_lock lock{ m_loggerMutex };
            m_logger = std::move(logger);
        }

        void setSampleRetentionCapacity(size_t capacity)
        {
            m_sampleRetentionCapacity = capacity;
        }

        std::optional<size_t> getSampleRetentionCapacity() const
        {
            return m_sampleRetentionCapacity;
        }

        void log(std::string message)
        {
            std::scoped_lock lock{ m_loggerMutex };
            m_logger(std::move(message));
        }

    private:
        static std::vector<InputSample> extractSamplesFromRetained(const std::optional<RingBuffer<InputSample>>& retainedSamples, double startTime, double endTime)
        {
            std::vector<InputSample> result{};
            if (retainedSamples.has_value())
            {
                for (const auto& sample : *retainedSamples)
                {
                    if (sample.Timestamp >= startTime && sample.Timestamp <= endTime)
                    {
                        result.push_back(sample);
                    }
                }
            }
            return result;
        }

        arcana::manual_dispatcher<256> m_callbackDispatcher{};
        std::optional<arcana::background_dispatcher<256>> m_processingDispatcher{};
        SchedulerT m_callbackScheduler{};
        SchedulerT m_processingScheduler{};
        std::vector<InputSample> m_samples{};
        std::vector<InputSample> m_processingSamples{};
        std::mutex m_samplesMutex{};
        TypedCollection m_dynamicSequenceProviders{};
        std::function<void(std::string)> m_logger{};
        std::mutex m_loggerMutex{};
        std::unordered_map<ContractId<>::IdT, std::unique_ptr<internal::CustomActionTypeOperations>> m_customActionIdsToOperations{};
        std::optional<size_t> m_sampleRetentionCapacity{};
    };
}
