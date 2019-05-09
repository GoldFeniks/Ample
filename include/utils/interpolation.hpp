#pragma once
#include <tuple>
#include <cstddef>
#include <utility>
#include <type_traits>
#include "types.hpp"
#include "utils.hpp"

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

        }// namespace __impl

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
                size_t i = 0;
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
                types::vector1d_t<T> res(n);
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

        template<typename T, typename V = T, typename I = linear_interpolation>
        class interpolated_data_1d {

        public:

            interpolated_data_1d() = delete;

            interpolated_data_1d(const interpolated_data_1d& other) {
                *this = other;
            }

            interpolated_data_1d(interpolated_data_1d&& other) {
                *this = std::move(other);
            }

            interpolated_data_1d& operator=(const interpolated_data_1d& other) {
                _interpolators = other._interpolators;
                _x = other._x;
                for (auto& it : _interpolators)
                    it._owner = this;
                return *this;
            }

            interpolated_data_1d& operator=(interpolated_data_1d&& other) {
                _interpolators = std::move(other._interpolators);
                _x = std::move(other._x);
                for (auto& it : _interpolators)
                    it._owner = this;
                return *this;
            }

            interpolated_data_1d(const types::vector1d_t<T>& x, const types::vector2d_t<V>& data) : _x(x) {
                _interpolators.reserve(data.size());
                for (const auto& it : data)
                    _interpolators.emplace_back(this, it);

            }
            interpolated_data_1d(const types::vector1d_t<T>& x, const types::vector1d_t<V>& data) : _x(x) {
                _interpolators.emplace_back(this, data);
            }

            interpolated_data_1d(types::vector1d_t<T>&& x, types::vector2d_t<V>&& data) : _x(std::move(x)) {
                _interpolators.reserve(data.size());
                for (const auto& it : data)
                    _interpolators.emplace_back(this, std::move(it));
            }
            interpolated_data_1d(types::vector1d_t<T>&& x, types::vector1d_t<V>&& data) : _x(std::move(x)) {
                _interpolators.emplace_back(this, std::move(data));
            }

            interpolated_data_1d(const types::vector1d_t<T>& x, types::vector2d_t<V>&& data) : _x(x) {
                _interpolators.reserve(data.size());
                for (const auto& it : data)
                    _interpolators.emplace_back(this, std::move(it));
            }
            interpolated_data_1d(const types::vector1d_t<T>& x, types::vector1d_t<V>&& data) : _x(x) {
                _interpolators.emplace_back(this, std::move(data));
            }

            const auto& operator[](const size_t i) const {
                return _interpolators[i];
            }

            auto size() const {
                return _interpolators.size();
            }

            const auto& x() const {
                return _x;
            }

        private:

            class interpolator {

            public:

                interpolator() = delete;

                interpolator(interpolated_data_1d* owner, const types::vector1d_t<V>& data) : _owner(owner), _data(data) {}

                interpolator(interpolated_data_1d* owner, types::vector1d_t<V>&& data) : _owner(owner), _data(std::move(data)) {}

                auto point(const T& x) const {
                    const auto [ix, jx] = utils::find_indices(_owner->_x, x);
                    return I::template line_point(_data[ix], _data[jx], _owner->_x[ix], _owner->_x[jx], x);
                }

                template<typename RV>
                auto line(const T& x0, const T& x1, RV& res) const {
                    I::template line(x0, x1, _owner->_x, _data, res);
                }

                template<typename RV>
                auto line(RV& res) const {
                    line(_owner->_x.front(), _owner->_x.back(), res);
                }

                auto line(const T& x0, const T& x1, const size_t n) const {
                    return I::template line(x0, x1, n, _owner->_x, _data);
                }

                auto line(const size_t n) const {
                    return line(_owner->_x.front(), _owner->_x.back(), n);
                }

                const auto& data() const {
                    return _data;
                }

                const auto& x() const {
                    return _owner->x();
                }

            private:

                interpolated_data_1d* _owner;
                const types::vector1d_t<V> _data;

                friend class interpolated_data_1d;

            };

            types::vector1d_t<T> _x;
            types::vector1d_t<interpolator> _interpolators;

        };

        template<typename T, typename V = T, typename I = linear_interpolation>
        class interpolated_data_2d {

        public:

            interpolated_data_2d() = delete;

            interpolated_data_2d(const interpolated_data_2d& other) {
                *this = other;
            }

            interpolated_data_2d(interpolated_data_2d&& other) {
                *this = std::move(other);
            }

            interpolated_data_2d& operator=(const interpolated_data_2d& other) {
                _interpolators = other._interpolators;
                _x = other._x;
                _y = other._y;
                for (auto& it : _interpolators)
                    it._owner = this;
                return *this;
            }

            interpolated_data_2d& operator=(interpolated_data_2d&& other) {
                _interpolators = std::move(other._interpolators);
                _x = std::move(other._x);
                _y = std::move(other._y);
                for (auto& it : _interpolators)
                    it._owner = this;
                return *this;
            }

            interpolated_data_2d(const types::vector1d_t<T>& x, const types::vector1d_t<T>& y,
                              const types::vector3d_t<V>& data) : _x(x), _y(y) {
                _interpolators.reserve(data.size());
                for (const auto& it : data)
                    _interpolators.emplace_back(this, it);

            }
            interpolated_data_2d(const types::vector1d_t<T>& x, const types::vector1d_t<T>& y,
                              const types::vector2d_t<V>& data) : _x(x), _y(y) {
                _interpolators.emplace_back(this, data);
            }

            interpolated_data_2d(types::vector1d_t<T>&& x, types::vector1d_t<T>&& y, types::vector3d_t<V>&& data) :
                _x(std::move(x)), _y(std::move(y)) {
                _interpolators.reserve(data.size());
                for (const auto& it : data)
                _interpolators.emplace_back(this, std::move(it));
            }
            interpolated_data_2d(types::vector1d_t<T>&& x, types::vector1d_t<T>&& y, types::vector2d_t<V>&& data) :
                    _x(std::move(x)), _y(std::move(y)) {
                _interpolators.emplace_back(this, std::move(data));
            }

            interpolated_data_2d(const types::vector1d_t<T>& x, const types::vector1d_t<T>& y, types::vector3d_t<V>&& data) :
                    _x(x), _y(y) {
                _interpolators.reserve(data.size());
                for (const auto& it : data)
                    _interpolators.emplace_back(this, std::move(it));
            }
            interpolated_data_2d(const types::vector1d_t<T>& x, const types::vector1d_t<T>& y, types::vector2d_t<V>&& data) :
                    _x(x), _y(y) {
                _interpolators.emplace_back(this, std::move(data));
            }

            const auto& operator[](const size_t i) const {
                return _interpolators[i];
            }

            auto size() const {
                return _interpolators.size();
            }

            const auto& x() const {
                return _x;
            }

            const auto& y() const {
                return _y;
            }

        private:

            class interpolator {

            public:

                interpolator() = delete;

                interpolator(interpolated_data_2d* owner, const types::vector2d_t<V>& data) : _owner(owner), _data(data) {}

                interpolator(interpolated_data_2d* owner, types::vector2d_t<V>&& data) : _owner(owner), _data(std::move(data)) {}

                auto point(const T& x, const T& y) const {
                    const auto [ix, jx] = utils::find_indices(_owner->_x, x);
                    const auto [iy, jy] = utils::find_indices(_owner->_y, y);
                    return I::template field_point(_data[ix][iy], _data[ix][jy], _data[jx][iy], _data[jx][jy],
                            _owner->_x[ix], _owner->_x[jx], _owner->_y[iy], _owner->_y[jy], x, y);
                }

                template<typename RV>
                auto line(const T& x, const T& y0, const T& y1, RV& res) const {
                    I::template field_line(x, y0, y1, _owner->_x, _owner->_y, _data, res);
                }

                template<typename RV>
                auto line(const T& x, RV& res) const {
                    line(x, _owner->_y.front(), _owner->_y.back(), res);
                }

                auto line(const T& x, const T& y0, const T& y1, const size_t n) const {
                    return I::template field_line(x, y0, y1, n, _owner->_x, _owner->_y, _data);
                }

                auto line(const T& x, const size_t n) const {
                    return line(x, _owner->_y.front(), _owner->_y.back(), n);
                }

                template<typename RV>
                auto field(const T& x0, const T& x1, const T& y0, const T& y1, RV& res) const {
                    I::template field(x0, x1, y0, y1, _owner->_x, _owner->_y, _data, res);
                }

                template<typename RV>
                auto field(RV& res) const {
                    return field(_owner->_x.front(), _owner->_x.back(), _owner->_y.front(), _owner->_y.back(), res);
                }

                auto field(const T& x0, const T& x1, const size_t nx, const T& y0, const T& y1, const size_t ny) const {
                    return I::template field(x0, x1, nx, y0, y1, ny, _owner->_x, _owner->_y, _data);
                }

                auto field(const size_t  nx, const size_t ny) const {
                    return field(_owner->_x.front(), _owner->_x.back(), nx, _owner->_y.front(), _owner->_y.back(), ny);
                }

                const auto& data() const {
                    return _data;
                }

                const auto& x() const {
                    return _owner->x();
                }

                const auto& y() const {
                    return _owner->y();
                }

            private:

                interpolated_data_2d* _owner;
                const types::vector2d_t<V> _data;

                friend class interpolated_data_2d;

            };

            types::vector1d_t<T> _x, _y;
            types::vector1d_t<interpolator> _interpolators;

        };


    }// namespace utils

}// namespace acstc
