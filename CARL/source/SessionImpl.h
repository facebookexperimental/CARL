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

namespace carl
{
    template<typename DescriptorT, size_t SequenceLength, typename = std::enable_if_t<std::is_trivially_copyable_v<DescriptorT>>>
    class DescriptorSequence
    {
        using WeakTableT = arcana::weak_table<typename Signal<gsl::span<const DescriptorT>>::HandlerT>;

    public:
        using DescriptorSignalT = Signal<gsl::span<const DescriptorT>>;

        class Provider : private WeakTableT, public DescriptorSignalT
        {
        public:
            Provider(Signal<gsl::span<const InputSample>>& signal)
                : DescriptorSignalT{ *static_cast<WeakTableT*>(this) }
                , m_ticket{ signal.addHandler([this](auto samples) { handleInputSample(samples); }) }
                {
                }

            void handleInputSample(gsl::span<const InputSample> samples)
            {
                for (size_t idx = 1; idx < samples.size(); ++idx)
                {
                    auto descriptor = DescriptorT::TryCreate(samples[idx], samples[idx - 1]);
                    if (descriptor.has_value())
                    {
                        m_sequence.emplace_back(std::move(descriptor.value()));
                    }
                }

                if (m_sequence.size() > SequenceLength)
                {
                    // Note: This does dangerous direct memory manipulation to avoid calling a bunch of
                    // range calls. This is allowable because of the requirement that descriptors be
                    // trivially copyable, but the behavior has the potential to be subtle and should
                    // be treated with caution.
                    m_buffer.resize(SequenceLength);
                    std::memcpy(m_buffer.data(), m_sequence.data() + (m_sequence.size() - SequenceLength), sizeof(DescriptorT) * SequenceLength);
                    m_buffer.swap(m_sequence);
                }

                auto& handlers = *static_cast<WeakTableT*>(this);
                handlers.apply_to_all([this](auto& callable) {
                    gsl::span<const DescriptorT> span{ m_sequence };
                    callable(span);
                });
            }

        private:
            const Signal<gsl::span<const InputSample>>::TicketT m_ticket;
            std::vector<DescriptorT> m_sequence{};
            std::vector<DescriptorT> m_buffer{};
        };
    };

    template <size_t SequenceLength, typename... DescriptorTs>
    class SessionImplBase
        : public arcana::weak_table<typename Signal<gsl::span<const InputSample>>::HandlerT>
        , public Signal<gsl::span<const InputSample>>
        , protected DescriptorSequence<DescriptorTs, SequenceLength>::Provider...
    {
    public:
        using SignalHandlersT = arcana::weak_table<typename Signal<gsl::span<const InputSample>>::HandlerT>;

        SessionImplBase()
            : Signal<gsl::span<const InputSample>>{ *static_cast<SignalHandlersT*>(this) }
            , DescriptorSequence<DescriptorTs, SequenceLength>::Provider{ *static_cast<Signal<gsl::span<const InputSample>>*>(this) }...
        {
        }

    protected:
        inline static constexpr size_t kSequenceLength{ SequenceLength };
    };

    class Session::Impl
        : public SessionImplBase<
        200,
        descriptor::HandPose<descriptor::Handedness::LeftHanded>,
        descriptor::HandGesture<descriptor::Handedness::LeftHanded>,
        descriptor::HandPose<descriptor::Handedness::RightHanded>,
        descriptor::HandGesture<descriptor::Handedness::RightHanded>,
        descriptor::TwoHandGesture>
    {
    public:
        Impl(double duration = 1. / 20.);

        static Session::Impl& getFromSession(Session& session);

        void addInputSample(const InputSample& inputSample);

        const double frameDuration{};

        template <typename DescriptorT>
        auto addHandler(std::function<void(gsl::span<const DescriptorT>)> handler)
        {
            return static_cast<typename DescriptorSequence<DescriptorT, kSequenceLength>::Provider*>(this)
                ->addHandler(std::move(handler));
        }

    private:
        std::vector<InputSample> m_samples{};
    };
}
