#pragma once
#include <cmath>
#include <tuple>
#include <cstddef>
#include <utility>
#include <functional>
#include <type_traits>
#include "types.hpp"
#include "utils.hpp"
#include "assert.hpp"
#include "delaunay.hpp"
#include "feniks/zip.hpp"

namespace ample::utils {

    namespace _impl {

        HAS_METHOD(prepare)

        template<typename T, typename C>
        auto fill_bilinear_interpolation_coefficients(const T& a, const T& b, const size_t n, const C& coords) {
            types::vector1d_t<size_t> is(n), js(n);
            types::vector1d_t<T> d1(n), d2(n), d12(n);
            const auto h = n > 1 ? (b - a) / (n - 1) : T(0);

            auto c = a;
            size_t i = 1, j = 0;

            while (c <= coords.front() && j < n) {
                d1[j] = is[j] = 0;
                d2[j] = d12[j] = 1;
                js[j] = (n > 1 ? 1 : 0);
                ++j;
                c += h;
            }

            while (c < coords.back() && j < n) {
                while (c > coords[i] && i < coords.size() - 1)
                    ++i;

                is[j] = i - 1;
                js[j] = i;

                d1[j] = c - coords[i - 1];
                d2[j] = coords[i] - c;
                d12[j] = coords[i] - coords[i - 1];

                ++j;
                c += h;
            }

            i = coords.size() - 1;
            while (j < n) {
                is[j] = i - 1;
                js[j] = i;

                d2[j] = 0;
                d1[j] = d12[j] = 1;
                ++j;
            }

            return std::make_tuple(is, js, d1, d2, d12);
        }

        template<typename T>
        struct remove_reference_wrapper {

            using type = T;

        };

        template<typename T>
        struct remove_reference_wrapper<std::reference_wrapper<T>> {

            using type = T;

        };

        template<typename T>
        using remove_reference_wrapper_t = typename remove_reference_wrapper<T>::type;

        template<typename T>
        auto& remove_reference_wrapper_v(T& value) {
            return static_cast<remove_reference_wrapper_t<T>&>(value);
        }

        class linear_interpolation {

        public:

            template<typename T, typename C>
            static auto line_point(const T& a, const T& b, const C& x0, const C& x1, const C& x) {
                return a + (b - a) * (x - x0) / (x1 - x0);
            }

            /**
                Interpolates value at point (x, y) inside a rectangle given values

                (x0, y0) -> a
                (x0, y1) -> b
                (x1, y0) -> c
                (x1, y1) -> d
            **/
            template<typename T, typename C>
            static auto field_point(const T& a, const T& b, const T& c, const T& d, const C& x0,
                                    const C& x1, const C& y0, const C& y1, const C& x, const C& y) {
                return line_point(line_point(a, b, y0, y1, y), line_point(c, d, y0, y1, y), x0, x1, x);
            }

            /**
                Interpolates value at point (x, y, z) inside a parallelogram given values

                (x0, y0, z0) -> a
                (x1, y0, z0) -> b
                (x0, y0, z1) -> c
                (x1, y0, z1) -> d
                (x0, y1, z0) -> e
                (x1, y1, z0) -> f
                (x0, y1, z1) -> g
                (x1, y1, z1) -> h

                  e------f
                 /|     /|
                a------b |
                | |    | |
                | g----|-h
                |/     |/
                c------d
            **/
            template<typename T, typename C>
            static auto area_point(
                    const T& a, const T& b, const T& c, const T& d, 
                    const T& e, const T& f, const T& g, const T& h,
                    const C& x0, const C& x1, 
                    const C& y0, const C& y1, 
                    const C& z0, const C& z1,
                    const C& x, const C& y, const C& z) {
                return line_point(
                    field_point(a, e, b, f, x0, x1, y0, y1, x, y),
                    field_point(c, g, d, h, x0, x1, y0, y1, x, y),
                    z0, z1, z
                );
            }

            template<typename T, typename C, typename V, typename It>
            static void line(const T& a, const T& b, const C& coords, const V& values, It begin, const It& end) {
                size_t i = 1;
                const auto n = std::distance(begin, end);
                const auto h = (b - a) / (n - 1);

                auto c = a;
                while (c <= coords.front() && begin != end) {
                    *(begin++) = values.front();
                    c += h;
                }

                while (c <= coords.back() && begin != end) {
                    while (c > coords[i] && i < coords.size() - 1)
                        ++i;

                    *(begin++) = values[i - 1] + (values[i] - values[i - 1]) * (c - coords[i - 1]) / (coords[i] - coords[i - 1]);
                    c += h;
                }

                while (begin != end)
                    *(begin++) = values.back();
            }

            template<typename T, typename C, typename V, typename RV>
            static void line(const T& a, const T& b, const C& coords, const V& values, RV& res) {
                line(a, b, coords, values, res.begin(), res.end());
            }

            template<typename T, typename C, typename V>
            static auto line(const T& a, const T& b, const size_t n, const C& coords, const V& values) {
                types::vector1d_t<std::decay_t<decltype(values[0])>> res(n);
                line(a, b, coords, values, res.begin(), res.end());
                return res;
            }

            template<typename T, typename C1, typename C2, typename V, typename RV>
            static auto field_line(const T& x, const T& y0, const T& y1, const C1& xs, const C2& ys, const V& values, RV& res) {
                const auto [xi, xj, dx1, dx2, dx12] = _impl::fill_bilinear_interpolation_coefficients(x, x, 1, xs);
                const auto [yi, yj, dy1, dy2, dy12] = _impl::fill_bilinear_interpolation_coefficients(y0, y1, res.size(), ys);
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
                const auto [xi, xj, dx1, dx2, dx12] = _impl::fill_bilinear_interpolation_coefficients(x0, x1, nx, xs);
                const auto [yi, yj, dy1, dy2, dy12] = _impl::fill_bilinear_interpolation_coefficients(y0, y1, ny, ys);
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

            template<typename T, typename C1, typename C2, typename C3, typename V, typename RV>
            static auto area_line(
                    const T& x, const T& y, 
                    const T& z0, const T& z1, 
                    const C1& xs, const C2& ys, const C3& zs, 
                    const V& values, RV& res) {
                types::vector2d_t<std::reference_wrapper<RV>> val{ types::vector1d_t<std::reference_wrapper<RV>>{ res } };
                area(x, x, y, y, z0, z1, xs, ys, zs, values, val);
            }

            template<typename T, typename C1, typename C2, typename C3, typename V, typename RV>
            static auto area_line(
                    const T& x, const T& y, 
                    const T& z0, const T& z1, const size_t& nz,
                    const C1& xs, const C2& ys, const C3& zs, 
                    const V& values) {
                types::vector1d_t<std::decay_t<decltype(values[0][0][0])>> res(nz);
                area_line(x, y, z0, z1, xs, ys, zs, values, res);
                return res;
            }

            template<typename T, typename C1, typename C2, typename C3, typename V, typename RV>
            static auto area_field(
                    const T& x, 
                    const T& y0, const T& y1,
                    const T& z0, const T& z1, 
                    const C1& xs, const C2& ys, const C3& zs, 
                    const V& values, RV& res) {
                types::vector1d_t<std::reference_wrapper<RV>> val{ res };
                area(x, x, y0, y1, z0, z1, xs, ys, zs, values, val);
            }

            template<typename T, typename C1, typename C2, typename C3, typename V, typename RV>
            static auto area_field(
                    const T& x, 
                    const T& y0, const T& y1, const size_t& ny, 
                    const T& z0, const T& z1, const size_t& nz,
                    const C1& xs, const C2& ys, const C3& zs, 
                    const V& values) {
                using value_type = std::decay_t<decltype(values[0][0][0])>;
                types::vector2d_t<value_type> res(ny, types::vector1d_t<value_type>(nz));
                area_field(x, y0, y1, z0, z1, xs, ys, zs, values, res);
                return res;
            }

            template<typename T, typename C1, typename C2, typename C3, typename V, typename RV>
            static auto area(
                    const T& x0, const T& x1,
                    const T& y0, const T& y1,
                    const T& z0, const T& z1, 
                    const C1& xs, const C2& ys, const C3& zs, 
                    const V& values, RV& res) {
                auto& result = remove_reference_wrapper_v(res);
                const auto nx = result.size();
                const auto ny = remove_reference_wrapper_v(result[0]).size();
                const auto nz = remove_reference_wrapper_v(remove_reference_wrapper_v(result[0])[0]).size();
                const auto [xi, xj, dx1, dx2, dx12] = _impl::fill_bilinear_interpolation_coefficients(x0, x1, nx, xs);
                const auto [yi, yj, dy1, dy2, dy12] = _impl::fill_bilinear_interpolation_coefficients(y0, y1, ny, ys);
                const auto [zi, zj, dz1, dz2, dz12] = _impl::fill_bilinear_interpolation_coefficients(z0, z1, nz, zs);
                for (size_t i = 0; i < nx; ++i) {
                    auto& ry = remove_reference_wrapper_v(result[i]);
                    for (size_t j = 0; j < ny; ++j) {
                        auto& rz = remove_reference_wrapper_v(ry[j]);
                        for (size_t k = 0; k < nz; ++k) {
                            rz[k] = 
                                (
                                    (
                                        (
                                            values[xi[i]][yi[j]][zi[k]] * dx2[i] +
                                            values[xj[i]][yi[j]][zi[k]] * dx1[i]
                                        ) * dy2[j]
                                        +
                                        (
                                            values[xi[i]][yj[j]][zi[k]] * dx2[i] +
                                            values[xj[i]][yj[j]][zi[k]] * dx1[i]
                                        ) * dy1[j]
                                    ) * dz2[k]
                                    +
                                    (
                                        (
                                            values[xi[i]][yi[j]][zj[k]] * dx2[i] +
                                            values[xj[i]][yi[j]][zj[k]] * dx1[i]
                                        ) * dy2[j]
                                        +
                                        (
                                            values[xi[i]][yj[j]][zj[k]] * dx2[i] +
                                            values[xj[i]][yj[j]][zj[k]] * dx1[i]
                                        ) * dy1[j]
                                    ) * dz1[k]
                                ) / dx12[i] / dy12[j] / dz12[k];
                        }
                    }
                }
            }

            template<typename T, typename C1, typename C2, typename C3, typename V>
            static auto area(
                    const T& x0, const T& x1, const size_t& nx, 
                    const T& y0, const T& y1, const size_t& ny,
                    const T& z0, const T& z1, const size_t& nz,
                    const C1& xs, const C2& ys, const C3& zs, 
                    const V& values) {
                using value_type = std::decay_t<decltype(values[0][0][0])>;
                types::vector3d_t<value_type> res(nx, types::vector2d_t<value_type>(ny, types::vector1d_t<value_type>(nz)));
                area(x0, x1, y0, y1, z0, z1, xs, ys, zs, values, res);
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
                    it._args = std::cref(_common_data);
                return *this;
            }

            interpolated_data& operator=(interpolated_data&& other) noexcept {
                std::swap(_interpolators, other._interpolators);
                std::swap(_common_data, other._common_data);
                for (auto& it : _interpolators)
                    it._args = std::cref(_common_data);
                return *this;
            }

            explicit interpolated_data(const Args&... args, const types::vector1d_t<data_t>& data) :
                _common_data(std::forward<const Args>(args)...) {
                prepare();
                _interpolators.reserve(data.size());
                for (const auto& it : data)
                    _interpolators.emplace_back(std::cref(_common_data), it);
            }

            explicit interpolated_data(const Args&... args, const data_t& data) :
                _common_data(std::forward<const Args>(args)...) {
                prepare();
                _interpolators.emplace_back(std::cref(_common_data), data);
            }

            explicit interpolated_data(const Args&... args, types::vector1d_t<data_t>&& data) :
                _common_data(std::forward<const Args>(args)...) {
                prepare();
                _interpolators.reserve(data.size());
                for (auto&& it : data)
                    _interpolators.emplace_back(std::cref(_common_data), std::move(it));
            }

            explicit interpolated_data(const Args&... args, data_t&& data) :
                _common_data(std::forward<const Args>(args)...) {
                prepare();
                _interpolators.emplace_back(std::cref(_common_data), std::move(data));
            }

            explicit interpolated_data(Args&&... args, types::vector1d_t<data_t>&& data) :
                _common_data(std::forward<Args>(args)...) {
                prepare();
                _interpolators.reserve(data.size());
                for (auto&& it : data)
                    _interpolators.emplace_back(std::cref(_common_data), std::move(it));
            }

            explicit interpolated_data(Args&&... args, data_t&& data) : _common_data(std::forward<Args>(args)...) {
                prepare();
                _interpolators.emplace_back(std::cref(_common_data), std::move(data));
            }

            auto& operator[](const size_t i) {
                return _interpolators[i];
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

            void erase_last(const size_t n = 0) {
                _interpolators.erase(_interpolators.end() - n, _interpolators.end());
            }

            void replace_data(const types::vector1d_t<data_t>& data) {
                utils::dynamic_assert(data.size() == _interpolators.size(),
                      "Incorrect number of data for interpolators. Expected ", _interpolators.size(), ", but got ", data.size());
                for (auto [interpolator, new_data] : feniks::zip(_interpolators, data))
                    interpolator.replace_data(new_data);
            }

            void replace_data(types::vector1d_t<data_t>&& data) {
                utils::dynamic_assert(data.size() == _interpolators.size(),
                                      "Incorrect number of data for interpolators. Expected ", _interpolators.size(), ", but got ", data.size());
                for (auto [interpolator, new_data] : feniks::zip(_interpolators, data))
                    interpolator.replace_data(std::move(new_data));
            }

        protected:

            std::tuple<Args...> _common_data;
            types::vector1d_t<I> _interpolators;

            void prepare() {
                if constexpr (has_prepare_v<I, std::tuple<Args...>>)
                    I::prepare(_common_data);
            }

        };

        template<typename I>
        class interpolated_data<I, std::tuple<>> {

        public:

            using data_t = typename I::data_t;

            interpolated_data() = default;

            explicit interpolated_data(const types::vector1d_t<data_t>& data) {
                _interpolators.reserve(data.size());
                for (const auto& it : data)
                    _interpolators.emplace_back(it);
            }

            explicit interpolated_data(const data_t& data) {
                _interpolators.emplace_back(data);
            }

            explicit interpolated_data(types::vector1d_t<data_t>&& data) {
                _interpolators.reserve(data.size());
                for (auto&& it : data)
                    _interpolators.emplace_back(std::move(it));
            }

            explicit interpolated_data(data_t&& data) {
                _interpolators.emplace_back(std::move(data));
            }

            const auto& operator[](const size_t i) const {
                return _interpolators[i];
            }

            void erase_last(const size_t n) {
                _interpolators.erase(_interpolators.end() - n, _interpolators.end());
            }

            auto size() const {
                return _interpolators.size();
            }

        protected:

            types::vector1d_t<I> _interpolators;

        };

        template<typename T, size_t N>
        auto get_vector_sizes(const types::vectornd_t<T, N>& vector) {
            const auto result = std::tuple<size_t>(vector.size());
            if constexpr (N == 1)
                return result;
            else
                return std::tuple_cat(result, get_vector_sizes<T, N - 1>(vector[0]));
        }

        template<size_t N, typename T, size_t K, typename S>
        void check_vector_sizes(const types::vectornd_t<T, K>& vector, const S& sizes) {
            utils::dynamic_assert(vector.size() == std::get<N - K>(sizes),
                  "Incorrect vector size at level ", N - K, ". Expected ", std::get<N - K>(sizes), ", but got ", vector.size());
            if constexpr (N > 1)
                for (const auto& it : vector)
                    check_vector_sizes<N, T, K - 1, S>(it, sizes);
        }

        template<typename T, size_t N>
        void check_vector_sizes(const types::vectornd_t<T, N>& a, const types::vectornd_t<T, N>& b) {
            const auto sizes = get_vector_sizes<T, N>(a);
            check_vector_sizes<N, T, N, decltype(sizes)>(b, sizes);
        }



    }// namespace _impl

    namespace interpolators {

        template<typename T, typename V>
        struct interpolator_1d {

            using line_t = types::vector1d_t<V>;

            virtual V point(const T& x) const = 0;
            virtual void line(const T& x0, const T& x1, line_t& res) const = 0;

            virtual line_t line(const T& x0, const T& x1, const size_t& n) const {
                line_t res(n);
                line(x0, x1, res);
                return res;
            }

        };

        template<typename T, typename V>
        struct interpolator_2d {

            using line_t = types::vector1d_t<V>;
            using field_t = types::vector2d_t<V>;

            virtual V point(const T& x, const T& y) const = 0;
            virtual void line(const T& x, const T& y0, const T& y1, line_t& res) const = 0;
            virtual void field(const T& x0, const T& x1, const T& y0, const T& y1, field_t& res) const = 0;

            virtual line_t line(const T& x, const T& y0, const T& y1, const size_t& n) const {
                line_t res(n);
                line(x, y0, y1, res);
                return res;
            }

            virtual field_t field(const T& x0, const T& x1, const size_t& nx, const T& y0, const T& y1, const size_t& ny) const {
                field_t res(nx, line_t(ny));
                field(x0, x1, y0, y1, res);
                return res;
            }

        };

        template<typename T, typename V>
        struct interpolator_3d {

            using line_t = types::vector1d_t<V>;
            using field_t = types::vector2d_t<V>;
            using area_t = types::vector3d_t<V>;

            virtual V point(const T& x, const T& y, const T& z) const = 0;
            virtual void line(const T& x, const T& y, const T& z0, const T& z1, line_t& res) const = 0;
            virtual void field(const T& x, const T& y0, const T& y1, const T& z0, const T& z1, field_t& res) const = 0;
            virtual void area(const T& x0, const T& x1, const T& y0, const T& y1, const T& z0, const T& z1, area_t& res) const = 0;

            virtual line_t line(const T& x, const T& y, const T& z0, const T& z1, const size_t& n) const {
                line_t res(n);
                line(x, y, z0, z1, res);
                return res;
            }

            virtual field_t field(const T& x, const T& y0, const T& y1, const size_t& ny, const T& z0, const T& z1, const size_t& nz) const {
                field_t res(ny, line_t(nz));
                field(x, y0, y1, z0, z1, res);
                return res;
            }

            virtual area_t area(
                    const T& x0, const T& x1, const size_t& nx,
                    const T& y0, const T& y1, const size_t& ny,
                    const T& z0, const T& z1, const size_t& nz) const {
                area_t res(nx, field_t(ny, line_t(nz)));
                area(x0, x1, y0, y1, z0, z1, res);
                return res;
            }

        };

        template<typename T, typename V = T>
        class linear_interpolator_1d : public interpolator_1d<T, V> {

        public:

            using data_t = types::vector1d_t<V>;
            using typename interpolator_1d<T, V>::line_t;
            using args_t = std::tuple<types::vector1d_t<T>>;

            using interpolator_1d<T, V>::line;

            linear_interpolator_1d(std::reference_wrapper<const args_t> args, const data_t& data) :
                    _args(std::move(args)), _data(data) {}
            linear_interpolator_1d(std::reference_wrapper<const args_t> args, data_t&& data) :
                    _args(std::move(args)), _data(std::move(data)) {}

            V point(const T& x) const override {
                const auto& xs = this->x();
                const auto [ix, jx] = find_indices(xs, x);
                return _impl::linear_interpolation::line_point(
                        _data[ix], _data[jx], xs[ix], xs[jx], x);
            }

            void line(const T& x0, const T& x1, line_t& res) const override {
                _impl::linear_interpolation::line(x0, x1, x(), _data, res);
            }

            void line(data_t& res) const {
                line(x()->front(), x()->back(), res);
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

            const auto& operator[](const size_t& i) const {
                return _data[i];
            }

            void replace_data(const data_t& data) {
                _impl::check_vector_sizes<V, 1>(_data, data);
                _data = data;
            }

            void replace_data(data_t&& data) {
                _impl::check_vector_sizes<V, 1>(_data, data);
                _data = std::move(data);
            }

        protected:

            data_t _data;
            std::reference_wrapper<const args_t> _args;

        private:

            template<typename, typename>
            friend class _impl::interpolated_data;

        };

        template<typename T, typename V = T>
        class linear_interpolator_2d : public interpolator_2d<T, V> {

        public:

            using data_t = types::vector2d_t<V>;
            using typename interpolator_2d<T, V>::line_t;
            using typename interpolator_2d<T, V>::field_t;
            using args_t = std::tuple<types::vector1d_t<T>, types::vector1d_t<T>>;

            using interpolator_2d<T, V>::line;
            using interpolator_2d<T, V>::field;

            linear_interpolator_2d(std::reference_wrapper<const args_t> args, const data_t& data) :
                    _args(std::move(args)), _data(data) {}
            linear_interpolator_2d(std::reference_wrapper<const args_t> args, data_t&& data) :
                    _args(std::move(args)), _data(std::move(data)) {}

            V point(const T& x, const T& y) const override {
                const auto& xs = this->x();
                const auto& ys = this->y();
                const auto [ix, jx] = find_indices(xs, x);
                const auto [iy, jy] = find_indices(ys, y);
                return _impl::linear_interpolation::field_point(
                        _data[ix][iy], _data[ix][jy], _data[jx][iy], _data[jx][jy],
                        xs[ix], xs[jx], ys[iy], ys[jy], x, y);
            }

            void line(const T& x, const T& y0, const T& y1, line_t& res) const override {
                _impl::linear_interpolation::field_line(x, y0, y1, this->x(), y(), _data, res);
            }

            void line(const T& x, line_t& res) const {
                line(x, y()->front(), y()->back(), res);
            }

            line_t line(const T& x, const size_t n) const {
                return line(x, y()->front(), y()->back(), n);
            }

            void field(const T& x0, const T& x1, const T& y0, const T& y1, field_t& res) const override {
                _impl::linear_interpolation::field(x0, x1, y0, y1, x(), y(), _data, res);
            }

            void field(data_t& res) const {
                field(x()->front(), x()->back(), y()->front(), y()->back(), res);
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

            const auto& operator[](const size_t& i) const {
                return _data[i];
            }

            void replace_data(const data_t& data) {
                _impl::check_vector_sizes<V, 2>(_data, data);
                _data = data;
            }

            void replace_data(data_t&& data) {
                _impl::check_vector_sizes<V, 2>(_data, data);
                _data = std::move(data);
            }

        protected:

            data_t _data;
            std::reference_wrapper<const args_t> _args;

        private:

            template<typename, typename>
            friend class _impl::interpolated_data;

        };

        template<typename T, typename V = T>
        class delaunay_interpolator_2d : public interpolator_2d<T, V> {

        public:

            using data_t = types::vector1d_t<V>;
            using typename interpolator_2d<T, V>::line_t;
            using typename interpolator_2d<T, V>::field_t;
            using args_t = std::tuple<delaunay_triangulation<T>>;

            using interpolator_2d<T, V>::line;
            using interpolator_2d<T, V>::field;

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

            void field(const T& x0, const T& x1, const T& y0, const T& y1, field_t& res) const {
                const auto hx = res.size() > 1 ? (x1 - x0) / (res.size() - 1) : T(0);
                for (size_t i = 0; i < res.size(); ++i)
                    line(x0 + i * hx, y0, y1, res[i]);
            }

            const auto& data() const {
                return _data;
            }

            void replace_data(const data_t& data) {
                _impl::check_vector_sizes<V, 1>(_data, data);
                _data = data;
            }

            void replace_data(data_t&& data) {
                _impl::check_vector_sizes<V, 1>(_data, data);
                _data = std::move(data);
            }

        private:

            data_t _data;
            std::reference_wrapper<const args_t> _args;

            template<typename C>
            V _point(const T& x, const T& y, const C& t) const {
                const auto& [a, b, c] = t.points();
                if (!t.contains_point({x, y})) {
                    const auto& [av, bv, cv] = t.get_points();
                    const auto da = std::hypot(av.x - x, av.y - y);
                    const auto db = std::hypot(bv.x - x, bv.y - y);
                    const auto dc = std::hypot(cv.x - x, cv.y - y);
                    if (da < db && da < dc)
                        return _data[a.index()];

                    if (db < dc)
                        return _data[b.index()];

                    return _data[c.index()];
                }

                const auto [u, v, w] = t.barycentric_coordinates({x, y});
                return u * _data[a.index()] + v * _data[b.index()] + w * _data[c.index()];
            }

            template<typename, typename>
            friend class _impl::interpolated_data;

        };

        template<typename T, typename V = T>
        class linear_interpolator_3d : public interpolator_3d<T, V> {

        public:

            using data_t = types::vector3d_t<V>;
            using typename interpolator_3d<T, V>::line_t;
            using typename interpolator_3d<T, V>::field_t;
            using typename interpolator_3d<T, V>::area_t;
            using args_t = std::tuple<types::vector1d_t<T>, types::vector1d_t<T>, types::vector1d_t<T>>;

            using interpolator_3d<T, V>::line;
            using interpolator_3d<T, V>::field;
            using interpolator_3d<T, V>::area;

            linear_interpolator_3d(std::reference_wrapper<const args_t> args, const data_t& data) :
                    _args(std::move(args)), _data(data) {}
            linear_interpolator_3d(std::reference_wrapper<const args_t> args, data_t&& data) :
                    _args(std::move(args)), _data(std::move(data)) {}

            V point(const T& x, const T& y, const T& z) const override {
                const auto& xs = this->x();
                const auto& ys = this->y();
                const auto& zs = this->z();
                const auto [ix, jx] = find_indices(xs, x);
                const auto [iy, jy] = find_indices(ys, y);
                const auto [iz, jz] = find_indices(zs, z);
                return _impl::linear_interpolation::area_point(
                        _data[ix][iy][iz], _data[jx][iy][iz], _data[ix][iy][jz], _data[jx][iy][jz],
                        _data[ix][jy][iz], _data[jx][jy][iz], _data[ix][jy][jz], _data[jx][jy][jz],
                        xs[ix], xs[jx], ys[iy], ys[jy], zs[iz], zs[jz],
                        x, y, z
                    );
            }

            void line(const T& x, const T& y, const T& z0, const T& z1, line_t& res) const override {
                _impl::linear_interpolation::area_line(x, y, z0, z1, this->x(), this->y(), this->z(), _data, res);
            }

            void line(const T& x, const T& y, line_t& res) const {
                line(x, y, z().front(), z().back(), res);
            }

            line_t line(const T& x, const T& y, const size_t n) const {
                return line(x, y, z().front(), z().back(), n);
            }

            void field(const T& x, const T& y0, const T& y1, const T& z0, const T& z1, field_t& res) const override {
                _impl::linear_interpolation::area_field(x, y0, y1, z0, z1, this->x(), y(), z(), _data, res);
            }

            void field(const T& x, field_t& res) const {
                field(x, y().front(), y().back(), z().front(), z().back(), res);
            }

            field_t field(const T& x, const size_t ny, const size_t nz) const {
                return interpolator_3d<T, V>::field(x, y().front(), y().back(), ny, z().front(), z().back(), nz);
            }

            void area(const T& x0, const T& x1, const T& y0, const T& y1, const T& z0, const T& z1, area_t& res) const override {
                _impl::linear_interpolation::area(x0, x1, y0, y1, z0, z1, x(), y(), z(), _data, res);
            }

            void area(area_t& res) const {
                area(x().front(), x().back(), y().front(), y().back(), z().front(), z().back(), res);
            }

            area_t area(const size_t& nx, const size_t& ny, const size_t& nz) const {
                return interpolator_3d<T, V>::area(
                    x().front(), x().back(), nx, 
                    y().front(), y().back(), ny,
                    z().front(), z().back(), nz);
            }

            inline const auto& x() const {
                return std::get<0>(_args.get());
            }

            inline const auto& y() const {
                return std::get<1>(_args.get());
            }

            inline const auto& z() const {
                return std::get<2>(_args.get());
            }

            inline const auto& data() const {
                return _data;
            }

            const auto& operator[](const size_t& i) const {
                return _data[i];
            }

            void replace_data(const data_t& data) {
                _impl::check_vector_sizes<V, 3>(_data, data);
                _data = data;
            }

            void replace_data(data_t&& data) {
                _impl::check_vector_sizes<V, 3>(_data, data);
                _data = std::move(data);
            }

        protected:

            data_t _data;
            std::reference_wrapper<const args_t> _args;

        private:

            template<typename, typename>
            friend class _impl::interpolated_data;

        };

    }// namespace interpolators

    template<typename I>
    using interpolated_data = _impl::interpolated_data<I, typename I::args_t>;

    template<typename T, typename V = T>
    using linear_interpolated_data_1d = interpolated_data<interpolators::linear_interpolator_1d<T, V>>;

    template<typename T, typename V = T>
    using linear_interpolated_data_2d = interpolated_data<interpolators::linear_interpolator_2d<T, V>>;
    
    template<typename T, typename V = T>
    using linear_interpolated_data_3d = interpolated_data<interpolators::linear_interpolator_3d<T, V>>;

    template<typename T, typename V = T>
    using delaunay_interpolated_data_2d = interpolated_data<interpolators::delaunay_interpolator_2d<T, V>>;

}// namespace ample::utils
