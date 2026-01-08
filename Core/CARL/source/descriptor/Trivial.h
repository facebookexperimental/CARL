/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "carl/Types.h"

#include <array>

namespace carl::descriptor
{
    // This namespace contains small, self-contained math constructs for use ONLY within
    // descriptors. These are needed because Eigen matrices are not trivially copyable,
    // so in order for descriptors to remain trivially copyable they need math primitives
    // that are also trivially copyable. Note that these types should be kept as minimal
    // as possible, containing only what's needed for descriptor comparison; Eigen should
    // still be used for all other math.
    namespace trivial
    {
        struct Point
        {
            NumberT X{};
            NumberT Y{};
            NumberT Z{};

            Point() = default;

            Point(const VectorT& v)
                : X{ v.x() }
                , Y{ v.y() }
                , Z{ v.z() }
            {
            }

            Point(NumberT x, NumberT y, NumberT z)
                : X{ x }
                , Y{ y }
                , Z{ z }
            {
            }

            Point& operator=(const VectorT& v)
            {
                X = v.x();
                Y = v.y();
                Z = v.z();
                return *this;
            }

            operator VectorT() const
            {
                return{ X, Y, Z };
            }

            Point operator-(const Point& other) const
            {
                return{ X - other.X, Y - other.Y, Z - other.Z };
            }

            NumberT lengthSq() const
            {
                return X * X + Y * Y + Z * Z;
            }

            NumberT length() const
            {
                return std::sqrt(lengthSq());
            }

            NumberT distanceSq(const Point& other) const
            {
                return (operator-(other)).lengthSq();
            }

            NumberT distance(const Point& other) const
            {
                return std::sqrt(distanceSq(other));
            }
        };

        struct Quaternion
        {
            NumberT X{};
            NumberT Y{};
            NumberT Z{};
            NumberT W{};

            Quaternion() = default;

            Quaternion(const QuaternionT& q)
                : X{ q.x() }
                , Y{ q.y() }
                , Z{ q.z() }
                , W{ q.w() }
            {
            }

            Quaternion(NumberT x, NumberT y, NumberT z, NumberT w)
                : X{ x }
                , Y{ y }
                , Z{ z }
                , W{ w }
            {
            }

            Quaternion& operator=(const QuaternionT& other)
            {
                X = other.x();
                Y = other.y();
                Z = other.z();
                W = other.w();
                return *this;
            }

            Quaternion operator*(const Quaternion& other) const
            {
                return{
                    W * other.X + X * other.W + Y * other.Z - Z * other.Y,
                    W * other.Y - X * other.Z + Y * other.W + Z * other.X,
                    W * other.Z + X * other.Y - Y * other.X + Z * other.W,
                    W * other.W - X * other.X - Y * other.Y - Z * other.Z,
                };
            }

            operator QuaternionT() const
            {
                return{ W, X, Y, Z };
            }

            Quaternion conjugate() const
            {
                return{ -X, -Y, -Z, W };
            }

            NumberT angularDistance(const Quaternion& other) const
            {
                auto qt{ operator*(other.conjugate()) };
                Point pt{ qt.X, qt.Y, qt.Z };
                return static_cast<NumberT>(2) * std::atan2(pt.length(), std::abs(qt.W));
            }
        };

        struct Transform
        {
            std::array<NumberT, 16> M{};

            Transform() = default;

            Transform(const TransformT& transform)
            {
                std::memcpy(M.data(), transform.matrix().data(), sizeof(M));
            }

            operator TransformT() const
            {
                auto m = M;
                return TransformT{ Eigen::Map<Eigen::Matrix4<NumberT>>{m.data()} };
            }

            Point operator *(const Point& p) const
            {
                return Point{
                    M[0] * p.X + M[4] * p.Y + M[8] * p.Z + M[12],
                    M[1] * p.X + M[5] * p.Y + M[9] * p.Z + M[13],
                    M[2] * p.X + M[6] * p.Y + M[10] * p.Z + M[14],
                };
            }
        };

        template<typename...>
        class Tuple;

        template<>
        class Tuple<> {};

        template<typename T, typename... Ts>
        class Tuple<T, Ts...>
        {
        public:
            Tuple() = default;

            Tuple(T&& t, Ts&&... ts)
                : m_t{ std::forward<T>(t) }, m_tail{ std::forward<Ts>(ts)... }
            {
            }

            template<size_t Idx, typename = std::enable_if_t<Idx != 0 && Idx <= sizeof...(Ts)>>
            auto& get()
            {
                return m_tail.template get<Idx - 1>();
            }

            template<size_t Idx, typename = std::enable_if_t<Idx != 0 && Idx <= sizeof...(Ts)>>
            const auto& get() const
            {
                return m_tail.template get<Idx - 1>();
            }

            template<size_t Idx, typename = std::enable_if_t<Idx == 0>>
            T& get()
            {
                return m_t;
            }

            template<size_t Idx, typename = std::enable_if_t<Idx == 0>>
            const T& get() const
            {
                return m_t;
            }

        private:
            T m_t;
            Tuple<Ts...> m_tail;
        };
    }
}
