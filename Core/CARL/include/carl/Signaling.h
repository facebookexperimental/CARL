/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <arcana/containers/weak_table.h>
#include <arcana/functional/inplace_function.h>

namespace carl
{
    /// <summary>
    /// Lightweight signaling primitive, allows consuming code to add handlers which will
    /// be called back when the signal fires.
    /// </summary>
    /// <typeparam name="...ArgsT">Arguments expected by handlers of this signal</typeparam>
    template <typename... ArgsT>
    class Signal
    {
        static constexpr size_t HANDLER_CAPACITY = 128;

    public:
        using HandlerT = stdext::inplace_function<void(ArgsT&...), HANDLER_CAPACITY, alignof(std::max_align_t), false>;
        using TicketT = typename arcana::weak_table<HandlerT>::ticket;

        Signal(arcana::weak_table<HandlerT>& handlers)
            : m_handlers{ handlers }
        {
        }

        Signal(const Signal&) = default;
        Signal(Signal&&) = default;
        Signal& operator=(const Signal&) = delete;
        Signal& operator=(Signal&&) = delete;

        template <typename CallableT>
        TicketT addHandler(CallableT&& callable)
        {
            return m_handlers.insert(std::forward<CallableT>(callable));
        }

    private:
        arcana::weak_table<HandlerT>& m_handlers;
    };
}
