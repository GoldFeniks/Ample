#pragma once
#include <tuple>
#include <cstddef>
#include "types.hpp"

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
                return line_point(y0, y1, line_point(x0, x1, a, b, x), line_point(x0, x1, c, d, x), y);
            }

            template<typename T, typename C, typename V>
            static auto line(const T& a, const T& b, const size_t n, const C& coords, const V& values) {
                size_t i = 0;
                const auto h = (b - a) / (n - 1);
                types::vector1d_t <T> res(n);
                for (size_t j = 0; j < n; ++j) {
                    const auto c = a + j * h;
                    while (c > coords[i] && i < coords.size() - 1)
                        ++i;
                    res[j] = values[i - 1] + (values[i] - values[i - 1]) * (c - coords[i - 1]) / (coords[i] - coords[i - 1]);
                }
                return res;
            }

            template<typename T, typename C1, typename C2, typename V>
            static auto field(const T& x0, const T& x1, const size_t nx, const T& y0, const T& y1, const size_t ny,
                                       const C1& xs, const C2& ys, const V& values) {
                const auto [xi, xj, dx1, dx2, dx12] = __impl::fill_bilinear_interpolation_coefficients(x0, x1, nx, xs);
                const auto [yi, yj, dy1, dy2, dy12] = __impl::fill_bilinear_interpolation_coefficients(y0, y1, ny, ys);
                types::vector2d_t <T> res(nx, types::vector1d_t<T>(ny));
                for (size_t i = 0; i < nx; ++i)
                    for (size_t j = 0; j < ny; ++j)
                        res[i][j] = (values[xi[i]][yi[j]] * dx2[i] * dy2[j] +
                                     values[xj[i]][yi[j]] * dx1[i] * dy2[j] +
                                     values[xi[i]][yj[j]] * dx2[i] * dy1[j] +
                                     values[xj[i]][yj[j]] * dx1[i] * dy1[j]) / dx12[i] / dy12[j];
                return res;
            }

        };

    }// namespace utils

}// namespace acstc
