/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "Types.h"

#include <Eigen/Geometry>

#include <gsl/span>

#include <array>
#include <cstring>
#include <optional>
#include <streambuf>
#include <string>
#include <vector>

namespace carl
{
    /// <summary>
    /// Helper class used to serialize data to byte buffers.
    /// </summary>
    class Serialization
    {
    public:
        Serialization(std::vector<uint8_t>& bytes)
            : m_bytes{ bytes }
        {
            m_bytes.clear();
        }

        template<typename T>
        void operator<<(const T& data)
        {
            size_t originalSize = m_bytes.size();
            m_bytes.resize(m_bytes.size() + sizeof(T));
            uint8_t* destination = m_bytes.data() + originalSize;
            std::memcpy(destination, &data, sizeof(T));
        }

        template<typename T>
        void operator<<(const std::vector<T>& data)
        {
            operator<<(static_cast<uint64_t>(data.size()));
            size_t originalSize = m_bytes.size();
            size_t vectorSizeBytes = data.size() * sizeof(T);
            m_bytes.resize(m_bytes.size() + vectorSizeBytes);
            uint8_t* destination = m_bytes.data() + originalSize;
            std::memcpy(destination, data.data(), vectorSizeBytes);
        }

        template<typename T, size_t Count>
        void operator<<(const std::array<T, Count>& data)
        {
            for (const auto& element : data)
            {
                operator<<(element);
            }
        }

        template<typename T>
        void operator<<(const std::optional<T>& data)
        {
            operator<<(data.has_value());
            if (data.has_value())
            {
                operator<<(data.value());
            }
        }

        template<typename T>
        static ptrdiff_t GetId(T& instance)
        {
            return reinterpret_cast<ptrdiff_t>(&instance);
        }

    private:
        std::vector<uint8_t>& m_bytes;
    };

    template<>
    inline void Serialization::operator<<<std::string>(const std::string& data)
    {
        size_t originalSize = m_bytes.size();
        m_bytes.resize(m_bytes.size() + data.size() + 1);
        uint8_t* destination = m_bytes.data() + originalSize;
        std::strcpy(reinterpret_cast<char*>(destination), data.c_str());
    }

    template<>
    inline void Serialization::operator<<<VectorT>(const VectorT& data)
    {
        operator<<(data.x());
        operator<<(data.y());
        operator<<(data.z());
    }

    template<>
    inline void Serialization::operator<<<QuaternionT>(const QuaternionT& data)
    {
        operator<<(data.x());
        operator<<(data.y());
        operator<<(data.z());
        operator<<(data.w());
    }

    template<>
    inline void Serialization::operator<<<TransformT>(const TransformT& data)
    {
        operator<<(QuaternionT{ data.rotation() });
        operator<<(VectorT{ data.translation() });
    }

    /// <summary>
    /// Helper class used to deserialize data from byte buffers.
    /// </summary>
    class Deserialization
    {
    public:
        Deserialization(const uint8_t* bytes, size_t length)
            : m_bytes{ bytes }
            , m_length{ length }
        {
        }

        Deserialization(gsl::span<const uint8_t> bytes)
            : m_bytes{ bytes.data() }
            , m_length{ bytes.size() }
        {
        }

        template<typename T>
        void operator>>(T& data)
        {
            assertAdequateLength(sizeof(T));

            data = *reinterpret_cast<const T*>(m_bytes);
            m_bytes += sizeof(T);
            m_length -= sizeof(T);
        }

        template<typename T>
        void operator>>(std::vector<T>& data)
        {
            uint64_t vectorSize64{};
            operator>>(vectorSize64);

            size_t vectorSize = static_cast<size_t>(vectorSize64);
            size_t vectorSizeBytes = vectorSize * sizeof(T);

            assertAdequateLength(vectorSizeBytes);

            data.resize(vectorSize);
            std::memcpy(data.data(), m_bytes, vectorSizeBytes);

            m_bytes += vectorSizeBytes;
            m_length -= vectorSizeBytes;
        }

        template<typename T, size_t Count>
        void operator>>(std::array<T, Count>& data)
        {
            for (size_t idx = 0; idx < Count; ++idx)
            {
                operator>>(data[idx]);
            }
        }

        template<typename T>
        void operator>>(std::optional<T>& data)
        {
            bool hasValue;
            operator>>(hasValue);
            if (hasValue)
            {
                data.emplace();
                operator>>(data.value());
            }
        }

        template<typename T>
        static T& GetReference(ptrdiff_t id)
        {
            return *reinterpret_cast<T*>(id);
        }

        size_t remainingBytes() const
        {
            return m_length;
        }

    private:
        const uint8_t* m_bytes{};
        size_t m_length{};

        void assertAdequateLength(size_t length)
        {
            if (length > m_length)
            {
                throw std::out_of_range{ "Deserialization failed: end of data reached" };
            }
        }
    };

    template<>
    inline void Deserialization::operator>><std::string>(std::string& data)
    {
        data = reinterpret_cast<const char*>(m_bytes);
        m_bytes += data.size() + 1;
        m_length -= data.size() + 1;
    }

    template<>
    inline void Deserialization::operator>><VectorT>(VectorT& data)
    {
        operator>>(data.x());
        operator>>(data.y());
        operator>>(data.z());
    }

    template<>
    inline void Deserialization::operator>><QuaternionT>(QuaternionT& data)
    {
        operator>>(data.x());
        operator>>(data.y());
        operator>>(data.z());
        operator>>(data.w());
    }

    template<>
    inline void Deserialization::operator>><TransformT>(TransformT& data)
    {
        QuaternionT rotation{};
        VectorT translation{};
        operator>>(rotation);
        operator>>(translation);
        data.fromPositionOrientationScale(translation, rotation, UNIT_SCALE);
    }
}
