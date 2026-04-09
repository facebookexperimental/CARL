/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <vector>
#include <cstddef>
#include <cassert>
#include <iterator>

namespace carl
{
    template<typename T>
    class RingBuffer
    {
    public:
        class Iterator
        {
        public:
            using iterator_category = std::random_access_iterator_tag;
            using value_type = T;
            using difference_type = std::ptrdiff_t;
            using pointer = const T*;
            using reference = const T&;

            Iterator() = default;

            reference operator*() const { return m_buffer->operator[](m_index); }
            pointer operator->() const { return &m_buffer->operator[](m_index); }
            reference operator[](difference_type n) const { return m_buffer->operator[](m_index + n); }

            Iterator& operator++() { ++m_index; return *this; }
            Iterator operator++(int) { auto tmp = *this; ++m_index; return tmp; }
            Iterator& operator--() { --m_index; return *this; }
            Iterator operator--(int) { auto tmp = *this; --m_index; return tmp; }

            Iterator& operator+=(difference_type n) { m_index += n; return *this; }
            Iterator& operator-=(difference_type n) { m_index -= n; return *this; }

            friend Iterator operator+(Iterator it, difference_type n) { it.m_index += n; return it; }
            friend Iterator operator+(difference_type n, Iterator it) { it.m_index += n; return it; }
            friend Iterator operator-(Iterator it, difference_type n) { it.m_index -= n; return it; }
            friend difference_type operator-(const Iterator& a, const Iterator& b) { return static_cast<difference_type>(a.m_index) - static_cast<difference_type>(b.m_index); }

            bool operator==(const Iterator& other) const { return m_index == other.m_index; }
            bool operator!=(const Iterator& other) const { return m_index != other.m_index; }
            bool operator<(const Iterator& other) const { return m_index < other.m_index; }
            bool operator<=(const Iterator& other) const { return m_index <= other.m_index; }
            bool operator>(const Iterator& other) const { return m_index > other.m_index; }
            bool operator>=(const Iterator& other) const { return m_index >= other.m_index; }

        private:
            friend class RingBuffer;
            Iterator(const RingBuffer* buffer, size_t index) : m_buffer{ buffer }, m_index{ index } {}

            const RingBuffer* m_buffer{};
            size_t m_index{};
        };

        explicit RingBuffer(size_t capacity)
            : m_data(capacity)
            , m_capacity{ capacity }
        {
            assert(capacity > 0);
        }

        void push_back(T value)
        {
            m_data[m_writePos] = std::move(value);
            m_writePos = (m_writePos + 1) % m_capacity;
            if (m_size < m_capacity)
            {
                ++m_size;
            }
        }

        const T& operator[](size_t idx) const
        {
            assert(idx < m_size);
            return m_data[(m_writePos + m_capacity - m_size + idx) % m_capacity];
        }

        T& operator[](size_t idx)
        {
            assert(idx < m_size);
            return m_data[(m_writePos + m_capacity - m_size + idx) % m_capacity];
        }

        const T& front() const { assert(m_size > 0); return operator[](0); }
        const T& back() const { assert(m_size > 0); return operator[](m_size - 1); }

        size_t size() const { return m_size; }
        size_t capacity() const { return m_capacity; }
        bool empty() const { return m_size == 0; }

        void clear()
        {
            m_size = 0;
            m_writePos = 0;
        }

        Iterator begin() const { return { this, 0 }; }
        Iterator end() const { return { this, m_size }; }

    private:
        std::vector<T> m_data{};
        size_t m_capacity{};
        size_t m_size{};
        size_t m_writePos{};
    };
}
