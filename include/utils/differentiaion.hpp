#pragma once
#include <algorithm>
#include <functional>
#include <type_traits>
#include "types.hpp"
#include "assert.hpp"
#include "feniks/zip.hpp"

namespace ample::utils {

    using namespace std::placeholders;

    namespace _impl {

        template<size_t N>
        struct x_derivative {

            template<typename T, typename F, typename... V>
            inline static void call(types::vectornd_t<T, N - 1>& v, F&& func, const V&... values) {
                for (size_t i = 0; i < v.size(); ++i)
                    x_derivative<N - 1>::template call<T, F, decltype(values[i])...>(v[i], std::forward<F>(func), values[i]...);
            }

        };

        template<>
        struct x_derivative<1> {

            template<typename T, typename F, typename... V>
            inline static void call(T& v, F&& func, const V&... values) {
                v = func(values...);
            }

        };

        template<typename T>
        struct step {

            using type = T;

        };

        template<typename T>
        struct step<std::complex<T>> {

            using type = T;

        };

        template<typename T>
        using step_t = typename step<T>::type;

        template<typename T, size_t M, size_t N>
        class differentiator_n_points_base {

        protected:

            using point_t = types::vectornd_t<T, M - 1>;

        public:

            static_assert(N >= 2, "Number of points for differentiation must be at least 2");

            void accept(point_t&& value) {
                if (_collected < N) {
                    _points[_collected] = std::move(value);
                }
                else {
                    _shift(std::make_index_sequence<N - 1>{});
                    _points.back() = std::move(value);
                }

                ++_collected;
            }

            void accept(const point_t& value) {
                point_t a(value);
                this->accept(std::move(a));
            }

        protected:

            size_t _collected = 0;
            std::array<point_t, N> _points{};

            template<size_t... I>
            void _shift(std::index_sequence<I...>) {
                (std::swap(_points[I], _points[I + 1]), ...);
            }

            [[nodiscard]] bool _ready() const {
                return _collected >= N;
            }

        };

    }// namespace _impl

    template<typename T, size_t N>
    class differentiator_3_points : public _impl::differentiator_n_points_base<T, N, 3> {

    private:

        using step_t = _impl::step_t<T>;
        using typename _impl::differentiator_n_points_base<T, N, 3>::point_t;

    public:

        static_assert(N >= 1 && N <= 3, "Number of dimensions must be between 1 and 3");

        explicit differentiator_3_points(std::array<step_t, N> d) : _d(std::move(d)) {}

        template<typename... V>
        explicit differentiator_3_points(const V&... values) : differentiator_3_points(std::array<step_t, N>{ static_cast<step_t>(values)... }) {
            static_assert(sizeof...(V) == N, "Incorrect number of arguments");
        }

        differentiator_3_points(std::initializer_list<step_t> values) {
            dynamic_assert(values.size() == N, "Incorrect number of arguments");
            for (auto [a, b] : feniks::zip(_d, values))
                a = b;
        }

        [[nodiscard]] bool next() {
            if (!_ready())
                return false;

            if (_returned < _collected) {
                ++_returned;
                return true;
            }

            return false;
        }

        const point_t& point() const {
            return _points[std::min(_returned - 1, size_t(2))];
        }

        const step_t& dx() const {
            return _d[0];
        }

        template<typename V = T>
        std::enable_if_t<(N > 1), const _impl::step_t<V>&> dy() const {
            return _d[1];
        }

        template<typename V = T>
        std::enable_if_t<(N > 2), const _impl::step_t<V>&> dz() const {
            return _d[2];
        }

        point_t x_derivative() const {
            dynamic_assert(_returned >= 1, "Not enough points");

            point_t result = _points[0];
            switch (_returned) {
                case 1:
                    _impl::x_derivative<N>::template call<T>(result, std::bind(_derivative_forward,  _1, _2, _3, dx()), _points[0], _points[1], _points[2]);
                    break;
                case 2:
                    _impl::x_derivative<N>::template call<T>(result, std::bind(_derivative_center,   _1, _2,     dx()), _points[0],             _points[2]);
                    break;
                default:
                    _impl::x_derivative<N>::template call<T>(result, std::bind(_derivative_backward, _1, _2, _3, dx()), _points[0], _points[1], _points[2]);
                    break;
            }

            return result;
        }

        template<typename V = T>
        std::enable_if_t<(N > 1), types::vectornd_t<V, N - 1>> y_derivative() const {
            dynamic_assert(_returned >= 1, "Not enough points");

            types::vectornd_t<T, N - 1> result = _points[0];
            const auto& points = point();
            const auto n = result.size();

            if constexpr (N == 2) {
                result[0] = _derivative_forward(points[0], points[1], points[2], dy());

                for (size_t i = 1; i < n - 1; ++i)
                    result[i] = _derivative_center(points[i - 1], points[i + 1], dy());

                result.back() = _derivative_backward(points[n - 3], points[n - 2], points[n - 1], dy());
            } else {
                for (size_t i = 0; i < result[0].size(); ++i)
                    result[0][i] = _derivative_forward(points[0][i], points[1][i], points[2][i], dy());

                for (size_t i = 1; i < n - 1; ++i)
                    for (size_t j = 0; j < result[i].size(); ++j)
                        result[i][j] = _derivative_center(points[i - 1][j], points[i + 1][j], dy());

                for (size_t i = 0; i < result.back().size(); ++i)
                    result.back()[i] = _derivative_backward(points[n - 3][i], points[n - 2][i], points[n - 1][i], dy());
            }

            return result;
        };

        template<typename V = T>
        std::enable_if_t<(N > 2), types::vectornd_t<V, N - 1>> z_derivative() const {
            dynamic_assert(_returned >= 1, "Not enough points");

            types::vectornd_t<T, N - 1> result = _points[0];
            const auto& points = point();

            for (size_t i = 0; i < result.size(); ++i) {
                const auto n = result[i].size();
                result[i][0] = _derivative_forward(points[i][0], points[i][1], points[i][2], dz());

                for (size_t j = 1; j < n - 1; ++j)
                    result[i][j] = _derivative_center(points[i][j - 1], points[i][j + 1], dz());

                result[i].back() = _derivative_backward(points[i][n - 3], points[i][n - 1], points[i][n - 1], dz());
            }

            return result;
        }

    private:

        using _impl::differentiator_n_points_base<T, N, 3>::_ready;
        using _impl::differentiator_n_points_base<T, N, 3>::_points;
        using _impl::differentiator_n_points_base<T, N, 3>::_collected;

        size_t _returned = 0;
        std::array<step_t, N> _d{};

        static T _derivative_forward(const T& a, const T& b, const T& c, const T& d) {
            return (-T(3) * a + T(4) * b - c) / (T(2) * d);
        }

        static T _derivative_center(const T& a, const T& c, const T& d) {
            return (c - a) / (T(2) * d);
        }

        static T _derivative_backward(const T& a, const T& b, const T& c, const T& d) {
            return (a - T(4) * b + T(3) * c) / (T(2) * d);
        }

    };

    template<typename T>
    using differentiator_3_points_1d = differentiator_3_points<T, 1>;

    template<typename T>
    using differentiator_3_points_2d = differentiator_3_points<T, 2>;

    template<typename T>
    using differentiator_3_points_3d = differentiator_3_points<T, 3>;

}// namespace ample::utils
