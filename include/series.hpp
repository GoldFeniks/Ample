#pragma once
#include <cmath>
#include <array>
#include <cstddef>
#include <type_traits>
#include "utils/types.hpp"

namespace ample::series {

    namespace _impl {

        template<typename T>
        struct coefficients {

            T a, b, c;

        };

        template<typename T>
        static constexpr auto pi = T(M_PI);

        template<typename T, size_t N>
        static constexpr auto den = T(2 * N) + T(1);

        template<typename T, size_t J, size_t N>
        static constexpr auto arg = J * pi<T> / den<T, N>;

        template<typename T, size_t J, size_t N>
        static const auto b = T(2) / den<T, N> * std::pow(std::sin(arg<T, J, N>), 2);

        template<typename T, size_t J, size_t N>
        static const auto c = std::pow(std::cos(arg<T, J, N>), 2);

        template<typename T, size_t N, size_t M, typename = std::enable_if_t<1 <= M>>
        struct setter {

            static void set(std::array<coefficients<T>, N>& coefficients) {
                setter<T, N, M - 1>::set(coefficients);
                coefficients[M - 1] = {T(0), b<T, M, N>, c<T, M, N>};
            }

        };

        template<typename T, size_t N>
        struct setter<T, N, 1> {

            static void set(std::array<coefficients<T>, N>& coefficients) {
                coefficients[0] = {T(1), b<T, 1, N> + c<T, 1, N>, c<T, 1, N>};
            }

        };

        template<typename T, size_t N>
        static auto get_coefficients() {
            std::array<coefficients<T>, N> result;
            setter<T, N, N>::set(result);
            return result;
        }

        template<size_t N, typename T = types::real_t>
        struct legendre_polynomial {

            template<typename V>
            static inline auto call(const V& value) {
                constexpr auto TN = T(N);
                return ((T(2) * TN - T(1)) * value * legendre_polynomial<N - 1, T>::call(value) -
                        (TN - T(1)) * legendre_polynomial<N - 2, T>::call(value)) / TN;
            }

        };

        template<typename T>
        struct legendre_polynomial<1, T> {

            template<typename V>
            static inline auto call(const V& value) {
                return value;
            }

        };

        template<typename T>
        struct legendre_polynomial<0, T> {

            template<typename V>
            static inline auto call(const V&) {
                return V(1);
            }

        };


    }// namespace _impl

    template<size_t N, typename T = types::real_t, typename = std::enable_if_t<0 != N>>
    class pade_series_coefficients {

    public:

        static constexpr auto size = N;

        pade_series_coefficients() : _coefficients(_impl::get_coefficients<T, N>()) {}

        explicit pade_series_coefficients(std::array<_impl::coefficients<T>, N> coefficients) :
            _coefficients(std::move(coefficients)) {}

        auto begin() {
            return _coefficients.begin();
        }

        auto begin() const {
            return _coefficients.begin();
        }

        auto end() {
            return _coefficients.end();
        }

        auto end() const {
            return _coefficients.end();
        }

        auto operator[](const size_t i) {
            return _coefficients[i];
        }

        auto operator[](const size_t i) const {
            return _coefficients[i];
        }

    private:

        const std::array<_impl::coefficients<T>, N> _coefficients;

    };

    template<size_t N, typename T = types::real_t>
    class legendre_polynomial {

    public:

        legendre_polynomial() = default;

        template<typename V>
        auto operator()(const V& value) const {
            return _impl::legendre_polynomial<N, T>::call(value);
        }

    };

}// namespace ample::series
