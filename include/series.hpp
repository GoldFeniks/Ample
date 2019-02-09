#pragma once
#include <cmath>
#include <cstddef>
#include "utils/types.hpp"

namespace acstc {

    namespace types = utils::types;

    namespace series {

        namespace impl {

            template<typename T>
            static constexpr auto pi = T(M_PI);

            template<typename T, size_t N>
            static constexpr auto den = T(2 * N) + T(1);

            template<typename T, size_t J, size_t N>
            static constexpr auto arg = J * pi<T> / den<T, N>;

            template<typename T, size_t J, size_t N>
            static constexpr auto a = T(2) / den<T, N> * std::pow(std::sin(arg<T, J, N>), 2);

            template<typename T, size_t J, size_t N>
            static constexpr auto b = std::pow(std::cos(arg<T, J, N>), 2);

            template<typename V, typename T, size_t N, size_t J = N>
            struct series {

                static auto call(const V& q) {
                    return a<T, J, N> * q / (T(1) + b<T, J, N> * q) + series<V, T, N, J - 1>::call(q);
                }

            };

            template<typename V, typename T, size_t N>
            struct series<V, T, N, 0> {

                static auto call(const V& q) {
                    return V(0);
                }

            };

        }// namespace impl

        template<size_t N, typename Val = types::complex_t, typename Arg = types::real_t>
        class pade_series {

        public:

            pade_series() = default;

            auto operator()(const Val& q) const {
                return Arg(1) + _series::call(q);
            }

        private:

            using _series = impl::series<Val, Arg, N>;

        };

    }// namespace series

}// namespace acstc