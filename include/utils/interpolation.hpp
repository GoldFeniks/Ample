#pragma once
#include <tuple>
#include <cstddef>
#include <utility>
#include <functional>
#include <type_traits>
#include "types.hpp"
#include "utils.hpp"
#include "delaunay.hpp"

namespace acstc {

    namespace utils {

        namespace __impl {

            template<typename T, typename C>
            auto fill_bilinear_interpolation_coefficients(const T& a, const T& b, const size_t n, const C& coords) {
                types::vector1d_t<size_t> is(n), js(n);
                types::vector1d_t<T> d1(n), d2(n), d12(n);
                const auto h = n > 1 ? (b - a) / (n - 1) : T(0);
                size_t i = 1;
                for (size_t j = 0; j < n; ++j) {
                    const auto c = a + j * h;
                    while (c > coords[i] && i < coords.size() - 1)
                        ++i;
                    is[j] = i - 1;
                    js[j] = i;
                    d1[j] = c - coords[i - 1];
                    d2[j] = coords[i] - c;
                    d12[j] = coords[i] - coords[i - 1];
                }
                return std::make_tuple(is, js, d1, d2, d12);
            }

            class linear_interpolation {

            public:

                template<typename T, typename C>
                static auto line_point(const T& a, const T& b, const C& x0, const C& x1, const C& x) {
                    return a + (b - a) * (x - x0) / (x1 - x0);
                }

                template<typename T, typename C>
                static auto field_point(const T& a, const T& b, const T& c, const T& d, const C& x0,
                                        const C& x1, const C& y0, const C& y1, const C& x, const C& y) {
                    return line_point(line_point(a, b, y0, y1, y), line_point(c, d, y0, y1, y), x0, x1, x);
                }

                template<typename T, typename C, typename V, typename RV>
                static auto line(const T& a, const T& b, const C& coords, const V& values, RV& res) {
                    size_t i = 1;
                    const auto h = (b - a) / (res.size() - 1);
                    for (size_t j = 0; j < res.size(); ++j) {
                        const auto c = a + j * h;
                        while (c > coords[i] && i < coords.size() - 1)
                            ++i;
                        res[j] = values[i - 1] + (values[i] - values[i - 1]) * (c - coords[i - 1]) / (coords[i] - coords[i - 1]);
                    }
                }

                template<typename T, typename C, typename V>
                static auto line(const T& a, const T& b, const size_t n, const C& coords, const V& values) {
                    types::vector1d_t<std::decay_t<decltype(values[0])>> res(n);
                    line(a, b, coords, values, res);
                    return res;
                }

                template<typename T, typename C1, typename C2, typename V, typename RV>
                static auto field_line(const T& x, const T& y0, const T& y1, const C1& xs, const C2& ys, const V& values, RV& res) {
                    const auto [xi, xj, dx1, dx2, dx12] = __impl::fill_bilinear_interpolation_coefficients(x, x, 1, xs);
                    const auto [yi, yj, dy1, dy2, dy12] = __impl::fill_bilinear_interpolation_coefficients(y0, y1, res.size(), ys);
                    for (size_t j = 0; j < res.size(); ++j)
                        res[j] = (values[xi[0]][yi[j]] * dx2[0] * dy2[j] +
                                  values[xj[0]][yi[j]] * dx1[0] * dy2[j] +
                                  values[xi[0]][yj[j]] * dx2[0] * dy1[j] +
                                  values[xj[0]][yj[j]] * dx1[0] * dy1[j]) / dx12[0] / dy12[j];
                }

                template<typename T, typename C1, typename C2, typename V>
                static auto field_line(const T& x, const T& y0, const T& y1, const size_t ny,
                                       const C1& xs, const C2& ys, const V& values) {
                    types::vector1d_t<std::decay_t<decltype(values[0][0])>> res(ny);
                    field_line(x, y0, y1, xs, ys, values, res);
                    return res;
                }

                template<typename T, typename C1, typename C2, typename V, typename RV>
                static auto field(const T& x0, const T& x1, const T& y0, const T& y1,
                                  const C1& xs, const C2& ys, const V& values, RV& res) {
                    const auto nx = res.size();
                    const auto ny = res[0].size();
                    const auto [xi, xj, dx1, dx2, dx12] = __impl::fill_bilinear_interpolation_coefficients(x0, x1, nx, xs);
                    const auto [yi, yj, dy1, dy2, dy12] = __impl::fill_bilinear_interpolation_coefficients(y0, y1, ny, ys);
                    for (size_t i = 0; i < nx; ++i)
                        for (size_t j = 0; j < ny; ++j)
                            res[i][j] = (values[xi[i]][yi[j]] * dx2[i] * dy2[j] +
                                         values[xj[i]][yi[j]] * dx1[i] * dy2[j] +
                                         values[xi[i]][yj[j]] * dx2[i] * dy1[j] +
                                         values[xj[i]][yj[j]] * dx1[i] * dy1[j]) / dx12[i] / dy12[j];
                }

                template<typename T, typename C1, typename C2, typename V>
                static auto field(const T& x0, const T& x1, const size_t nx, const T& y0, const T& y1, const size_t ny,
                                  const C1& xs, const C2& ys, const V& values) {
                    types::vector2d_t<std::decay_t<decltype(values[0][0])>>
                            res(nx, types::vector1d_t<std::decay_t<decltype(values[0][0])>>(ny));
                    field(x0, x1, y0, y1, xs, ys, values, res);
                    return res;
                }

            };

            template<typename I, typename T>
            class interpolated_data;

            template<typename I, typename... Args>
            class interpolated_data<I, std::tuple<Args...>> {

            public:

                using data_t = typename I::data_t;

                interpolated_data() = default;

                interpolated_data(const interpolated_data& other) {
                    *this = other;
                }

                interpolated_data(interpolated_data&& other) noexcept {
                    *this = std::move(other);
                }

                interpolated_data& operator=(const interpolated_data& other) {
                    _interpolators = other._interpolators;
                    _common_data = other._common_data;
                    for (auto& it : _interpolators)
                        it.data_ref = std::cref(_common_data);
                    return *this;
                }

                interpolated_data& operator=(interpolated_data&& other) noexcept {
                    std::swap(_interpolators, other._interpolators);
                    std::swap(_common_data, other._common_data);
                    for (auto& it : _interpolators)
                        it._args = std::cref(_common_data);
                    return *this;
                }

                interpolated_data(const Args&... args, const types::vector1d_t<data_t>& data) :
                    _common_data(std::forward<const Args>(args)...) {
                    _interpolators.reserve(data.size());
                    for (const auto& it : data)
                        _interpolators.emplace_back(std::cref(_common_data), it);
                }

                interpolated_data(const Args&... args, const data_t& data) :
                    _common_data(std::forward<const Args>(args)...) {
                    _interpolators.emplace_back(std::cref(_common_data), data);
                }

                interpolated_data(const Args&... args, types::vector1d_t<data_t>&& data) :
                    _common_data(std::forward<const Args>(args)...) {
                    _interpolators.reserve(data.size());
                    for (auto&& it : data)
                        _interpolators.emplace_back(std::cref(_common_data), std::move(it));
                }

                interpolated_data(const Args&... args, data_t&& data) :
                    _common_data(std::forward<const Args>(args)...) {
                    _interpolators.emplace_back(std::cref(_common_data), std::move(data));
                }

                interpolated_data(Args&&... args, types::vector1d_t<data_t>&& data) :
                    _common_data(std::forward<Args>(args)...) {
                    _interpolators.reserve(data.size());
                    for (auto&& it : data)
                        _interpolators.emplace_back(std::cref(_common_data), std::move(it));
                }

                interpolated_data(Args&&... args, data_t&& data) : _common_data(std::forward<Args>(args)...) {
                    _interpolators.emplace_back(std::cref(_common_data), std::move(data));
                }

                const auto& operator[](const size_t i) const {
                    return _interpolators[i];
                }

                auto size() const {
                    return _interpolators.size();
                }

                template<size_t N>
                const auto& get() const {
                    return std::get<N>(_common_data);
                }

            protected:

                std::tuple<Args...> _common_data;
                types::vector1d_t<I> _interpolators;

            };

            template<typename I>
            class interpolated_data<I, std::tuple<>> {

            public:

                using data_t = typename I::data_t;

                interpolated_data() = default;

                interpolated_data(const types::vector1d_t<data_t>& data) {
                    _interpolators.reserve(data.size());
                    for (const auto& it : data)
                        _interpolators.emplace_back(it);
                }

                interpolated_data(const data_t& data) {
                    _interpolators.emplace_back(data);
                }

                interpolated_data(types::vector1d_t<data_t>&& data) {
                    _interpolators.reserve(data.size());
                    for (auto&& it : data)
                        _interpolators.emplace_back(std::move(it));
                }

                interpolated_data(data_t&& data) {
                    _interpolators.emplace_back(std::move(data));
                }

                const auto& operator[](const size_t i) const {
                    return _interpolators[i];
                }

                auto size() const {
                    return _interpolators.size();
                }

            protected:

                types::vector1d_t<I> _interpolators;

            };

        }// namespace __impl

        namespace interpolators {

            template<typename T, typename V>
            struct interpolator_1d {

                using line_t = types::vector1d_t<V>;

                virtual V point(const T& x) const = 0;
                virtual void line(const T& x0, const T& x1, line_t& res) const = 0;

            };

            template<typename T, typename V>
            struct interpolator_2d {

                using line_t = types::vector1d_t<V>;
                using field_t = types::vector2d_t<V>;

                virtual V point(const T& x, const T& y) const = 0;
                virtual void line(const T& x, const T& y0, const T& y1, line_t& res) const = 0;
                virtual void field(const T& x0, const T& x1, const T& y0, const T& y1, field_t& res) const = 0;

            };

            template<typename T, typename V = T>
            class linear_interpolator_1d : interpolator_1d<T, V> {

            public:

                using data_t = types::vector1d_t<V>;
                using typename interpolator_1d<T, V>::line_t;
                using args_t = std::tuple<types::vector1d_t<T>>;

                linear_interpolator_1d(std::reference_wrapper<const args_t> args, const data_t& data) :
                        _args(std::move(args)), _data(data) {}
                linear_interpolator_1d(std::reference_wrapper<const args_t> args, data_t&& data) :
                        _args(std::move(args)), _data(std::move(data)) {}

                V point(const T& x) const override {
                    const auto& xs = this->x();
                    const auto [ix, jx] = utils::find_indices(xs, x);
                    return __impl::linear_interpolation::line_point(
                            _data[ix], _data[jx], xs[ix], xs[jx], x);
                }

                void line(const T& x0, const T& x1, line_t& res) const override {
                    __impl::linear_interpolation::line(x0, x1, x(), _data, res);
                }

                void line(data_t& res) const {
                    line(x()->front(), x()->back(), res);
                }

                auto line(const T& x0, const T& x1, const size_t n) const {
                    line_t res(n);
                    line(x0, x1, res);
                    return res;
                }

                auto line(const size_t n) const {
                    return line(x().front(), x().back(), n);
                }

                inline const auto& x() const {
                    return std::get<0>(_args.get());
                }

                inline const auto& data() const {
                    return _data;
                }

            protected:

                data_t _data;
                std::reference_wrapper<const args_t> _args;

            private:

                template<typename, typename>
                friend class __impl::interpolated_data;

            };

            template<typename T, typename V = T>
            class linear_interpolator_2d : interpolator_2d<T, V> {

            public:

                using data_t = types::vector2d_t<V>;
                using typename interpolator_2d<T, V>::line_t;
                using typename interpolator_2d<T, V>::field_t;
                using args_t = std::tuple<types::vector1d_t<T>, types::vector1d_t<T>>;

                linear_interpolator_2d(std::reference_wrapper<const args_t> args, const data_t& data) :
                        _args(std::move(args)), _data(data) {}
                linear_interpolator_2d(std::reference_wrapper<const args_t> args, data_t&& data) :
                        _args(std::move(args)), _data(std::move(data)) {}

                V point(const T& x, const T& y) const override {
                    const auto& xs = this->x();
                    const auto& ys = this->y();
                    const auto [ix, jx] = utils::find_indices(xs, x);
                    const auto [iy, jy] = utils::find_indices(ys, y);
                    return __impl::linear_interpolation::field_point(
                            _data[ix][iy], _data[ix][jy], _data[jx][iy], _data[jx][jy],
                            xs[ix], xs[jx], ys[iy], ys[jy], x, y);
                }

                void line(const T& x, const T& y0, const T& y1, line_t& res) const override {
                    __impl::linear_interpolation::field_line(x, y0, y1, this->x(), y(), _data, res);
                }

                void line(const T& x, line_t& res) const {
                    line(x, y()->front(), y()->back(), res);
                }

                line_t line(const T& x, const T& y0, const T& y1, size_t n) const {
                    line_t res(n);
                    line(x, y0, y1, res);
                    return res;
                }

                line_t line(const T& x, const size_t n) const {
                    return line(x, y()->front(), y()->back(), n);
                }

                void field(const T& x0, const T& x1, const T& y0, const T& y1, field_t& res) const override {
                    __impl::linear_interpolation::field(x0, x1, y0, y1, x(), y(), _data, res);
                }

                void field(data_t& res) const {
                    field(x()->front(), x()->back(), y()->front(), y()->back(), res);
                }

                field_t field(const T& x0, const T& x1, const size_t nx, const T& y0, const T& y1, const size_t ny) const {
                    field_t res(nx, line_t(ny));
                    field(x0, x1, y0, y1, res);
                    return res;
                }

                field_t field(const size_t  nx, const size_t ny) const {
                    return field(x().front(), x().back(), nx, y().front(), y().back(), ny);
                }

                inline const auto& x() const {
                    return std::get<0>(_args.get());
                }

                inline const auto& y() const {
                    return std::get<1>(_args.get());
                }

                inline const auto& data() const {
                    return _data;
                }

            protected:

                data_t _data;
                std::reference_wrapper<const args_t> _args;

            private:

                template<typename, typename>
                friend class __impl::interpolated_data;

            };

            template<typename T, typename V = T>
            class delaunay_interpolator_2d : interpolator_2d<T, V> {

            public:

                using data_t = types::vector1d_t<V>;
                using typename interpolator_2d<T, V>::line_t;
                using typename interpolator_2d<T, V>::field_t;
                using args_t = std::tuple<delaunay_triangulation<T>>;

                delaunay_interpolator_2d(std::reference_wrapper<const args_t> args, const data_t& data) :
                        _args(std::move(args)), _data(data) {}

                delaunay_interpolator_2d(std::reference_wrapper<const args_t> args, data_t&& data) :
                        _args(std::move(args)), _data(std::move(data)) {}

                V point(const T& x, const T& y) const override {
                    return _point(x, y, std::get<0>(_args.get()).find_triangle(x, y).get());
                }

                void line(const T& x, const T& y0, const T& y1, line_t& res) const override {
                    const auto h = res.size() > 1 ? (y1 - y0) / (res.size() - 1) : T(0);
                    auto t = std::get<0>(_args.get()).find_triangle(x, y0);
                    res[0] = _point(x, y0, t.get());
                    for (size_t i = 1; i < res.size(); ++i) {
                        t = std::get<0>(_args.get()).find_triangle(t, x, y0 + h * i);
                        res[i] = _point(x, y0 + h * i, t.get());
                    }
                }

                line_t line(const T& x, const T& y0, const T& y1, const size_t nx) const {
                    line_t res(nx);
                    line(x, y0, y1, res);
                    return res;
                }

                void field(const T& x0, const T& x1, const T& y0, const T& y1, field_t& res) const {
                    const auto hx = res.size() > 1 ? (x1 - x0) / (res.size() - 1) : T(0);
                    for (size_t i = 0; i < res.size(); ++i)
                        line(x0 + i * hx, y0, y1, res[i]);
                }

                field_t field(const T& x0, const T& x1, const size_t nx, const T& y0, const T& y1, const size_t ny) const {
                    field_t res(nx, line_t(ny));
                    field(x0, x1, y0, y1, res);
                    return res;
                }

            private:

                data_t _data;
                std::reference_wrapper<const args_t> _args;

                template<typename C>
                V _point(const T& x, const T& y, const C& t) const {
                    const auto& [a, b, c] = t.points();
                    const auto& [o, p, q] = t.get_points();
                    const auto  [u, v, w] = t.barycentric_coordinates({x, y});
                    const auto buff = u * _data[a.index()] + v * _data[b.index()] + w * _data[c.index()];
                    return u * _data[a.index()] + v * _data[b.index()] + w * _data[c.index()];
                }

                template<typename, typename>
                friend class __impl::interpolated_data;

            };

        }// namespace interpolators

        template<typename I>
        using interpolated_data = __impl::interpolated_data<I, typename I::args_t>;

        template<typename T, typename V = T>
        using linear_interpolated_data_1d = interpolated_data<interpolators::linear_interpolator_1d<T, V>>;

        template<typename T, typename V = T>
        using linear_interpolated_data_2d = interpolated_data<interpolators::linear_interpolator_2d<T, V>>;

        template<typename T, typename V = T>
        using delaunay_interpolated_data_2d = interpolated_data<interpolators::delaunay_interpolator_2d<T, V>>;

    }// namespace utils

}// namespace acstc
