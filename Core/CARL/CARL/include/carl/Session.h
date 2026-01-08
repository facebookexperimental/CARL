/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "carl/ContractId.h"
#include "carl/InputSample.h"
#include "carl/Example.h"

#include <arcana/threading/dispatcher.h>

#include <memory>

namespace carl
{
    namespace internal
    {
        struct CustomActionTypeOperations
        {
            using T = std::shared_ptr<void>;

            using TryCreateFunctionT = std::function<std::optional<T>(const InputSample& current, const InputSample& prior)>;
            TryCreateFunctionT TryCreate{};

            using DistanceFunctionT = std::function<NumberT(const T& a, const T& a0, const T& b, const T& b0, gsl::span<const NumberT>)>;
            DistanceFunctionT Distance{};

            using LerpFunctionT = std::function<T(const T& a, const T& b, NumberT t)>;
            LerpFunctionT Lerp{};

            using CalculateTuningFunctionT = std::function<std::array<NumberT, 32>(gsl::span<const action::Example> examples)>;
            CalculateTuningFunctionT CalculateTuning{};


            template<typename T>
            static std::shared_ptr<void> toOwningTypeErasedPtr(T data)
            {
                constexpr auto deleter = [](void* p) { delete reinterpret_cast<T*>(p); };
                return{ reinterpret_cast<void*>(new T(std::move(data))), deleter };
            }

            template<typename T>
            static const T& fromOwningTypeErasedPtr(const std::shared_ptr<void>& ptr)
            {
                return *reinterpret_cast<const T*>(ptr.get());
            }
        };
    }

    template<typename T>
    struct TemplatedCustomActionTypeOperations
    {
        using TryCreateFunctionT = std::function<std::optional<T>(const InputSample& current, const InputSample& prior)>;
        TryCreateFunctionT TryCreate{};

        using DistanceFunctionT = std::function<NumberT(const T& a, const T& a0, const T& b, const T& b0, gsl::span<const NumberT>)>;
        DistanceFunctionT Distance{};

        using LerpFunctionT = std::function<T(const T& a, const T& b, NumberT t)>;
        LerpFunctionT Lerp{};

        using CalculateTuningFunctionT = std::function<std::array<NumberT, 32>(gsl::span<const action::Example> examples)>;
        CalculateTuningFunctionT CalculateTuning{};
    };

    /// <summary>
    /// CARL's characterization of a timespan and usage in which actions can be recognized.
    /// </summary>
    class Session
    {
    public:
        class Impl;
        using SchedulerT = stdext::inplace_function<void(stdext::inplace_function<void(), 128>&&), 128>;

        Session(bool singleThreaded = false);
        ~Session();

        void addInput(const InputSample&);

        SchedulerT& callbackScheduler();
        SchedulerT& processingScheduler();
        void setLogger(std::function<void(std::string)> logger);
        void log(std::string message);
        void tickCallbacks(arcana::cancellation& token);

        template<typename T>
        ContractId<>::IdT enableCustomActionType(TemplatedCustomActionTypeOperations<T> operations)
        {
            using InternalT = internal::CustomActionTypeOperations::T;

            auto opsPtr = std::make_shared<decltype(operations)>(std::move(operations));

            internal::CustomActionTypeOperations internalOps{};

            internalOps.TryCreate = [opsPtr](const InputSample& current, const InputSample& prior) -> std::optional<InternalT> {
                auto created = opsPtr->TryCreate(current, prior);
                if (created.has_value())
                {
                    return internal::CustomActionTypeOperations::toOwningTypeErasedPtr(std::move(created.value()));
                }
                else
                {
                    return{};
                }
                };

            internalOps.Distance = [opsPtr](const InternalT& a, const InternalT& a0, const InternalT& b, const InternalT& b0, gsl::span<const NumberT> tuning) -> NumberT {
                const auto& aRef = internal::CustomActionTypeOperations::fromOwningTypeErasedPtr<T>(a);
                const auto& a0Ref = internal::CustomActionTypeOperations::fromOwningTypeErasedPtr<T>(a0);
                const auto& bRef = internal::CustomActionTypeOperations::fromOwningTypeErasedPtr<T>(b);
                const auto& b0Ref = internal::CustomActionTypeOperations::fromOwningTypeErasedPtr<T>(b0);
                return opsPtr->Distance(aRef, a0Ref, bRef, b0Ref, tuning);
                };

            internalOps.Lerp = [opsPtr](const InternalT& a, const InternalT& b, NumberT t) -> InternalT {
                const auto& aRef = internal::CustomActionTypeOperations::fromOwningTypeErasedPtr<T>(a);
                const auto& bRef = internal::CustomActionTypeOperations::fromOwningTypeErasedPtr<T>(b);
                return internal::CustomActionTypeOperations::toOwningTypeErasedPtr(opsPtr->Lerp(aRef, bRef, t));
                };

            internalOps.CalculateTuning = opsPtr->CalculateTuning;

            return internalEnableCustomActionType(std::move(internalOps));
        }

    private:
        friend class Impl;
        std::unique_ptr<Impl> m_impl{};

        ContractId<>::IdT internalEnableCustomActionType(internal::CustomActionTypeOperations);
    };
}
