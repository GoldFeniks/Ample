#pragma once

#include <cmath>
#include <cstddef>
#include "dork.hpp"
#include "utils/types.hpp"
#include "utils/interpolation.hpp"

namespace acstc {

    namespace rays {

        namespace __impl {

            template<typename T, typename V>
            V calc_derivative_forward(const V& y0, const V& y1, const V& y2, const T& x0, const T& x1, const T& x2) {
                const T h = x1 - x0, d = x2 - x0;
                return (d * d * (y1 - y0) + h * h * (y0 - y2)) / (d * (x2 - x1) * h);
            }

            template<typename T, typename V>
            V calc_derivative_backward(const V& y0, const V& y1, const V& y2, const T& x0, const T& x1, const T& x2) {
                const T h = x2 - x1, d = x2 - x0;
                return (d * d * (y2 - y1) + h * h * (y0 - y2)) / (d * (x1 - x0) * h);
            }

            template<typename T, typename V>
            auto calc_derivative(const utils::linear_interpolated_data_1d<T, V>& values) {
                const auto& x = values.template get<0>();
                types::vector2d_t<V> result(values.size());
                for (size_t j = 0; j < values.size(); ++j) {
                    const auto& data = values[j].data();
                    result[j].resize(x.size());
                    for (size_t i = 1; i < x.size() - 1; ++i)
                            result[j][i] = (data[i + 1] - data[i - 1]) / (x[i + 1] - x[i - 1]);
                    result[j][0] = calc_derivative_forward(data[0], data[1], data[2], x[0], x[1], x[2]);
                    result[j].back() = calc_derivative_backward(data[x.size() - 3], data[x.size() - 2], data.back(),
                            x[x.size() - 3], x[x.size() - 2], x.back());
                }
                return utils::linear_interpolated_data_1d<T, V>(x, result);
            }

            template<typename T, typename V>
            auto calc_derivatives(const acstc::utils::linear_interpolated_data_2d<T, V>& values) {
                const auto& x = values.template get<0>();
                const auto& y = values.template get<1>();
                types::vector3d_t<V> resx(values.size()), resy(values.size());
                for (size_t j = 0; j < values.size(); ++j) {
                    const auto& data = values[j].data();
                    resx[j].resize(x.size(), types::vector1d_t<V>(y.size()));
                    resy[j].resize(x.size(), types::vector1d_t<V>(y.size()));
                    for (size_t i = 1; i < x.size() - 1; ++i) {
                            resx[j][i][0] = (data[i + 1][0] - data[i - 1][0]) / (x[i + 1] - x[i - 1]);
                            resx[j][i].back() = (data[i + 1].back() - data[i - 1].back()) / (x[i + 1] - x[i - 1]);
                            for (size_t k = 1; k < y.size() - 1; ++k) {
                                    resx[j][i][k] = (data[i + 1][k] - data[i - 1][k]) / (x[i + 1] - x[i - 1]);
                                    resy[j][i][k] = (data[i][k + 1] - data[i][k - 1]) / (y[k + 1] - y[k - 1]);
                            }
                    }
                    for (size_t i = 0; i < x.size(); ++i) {
                            resy[j][i][0] = calc_derivative_forward(data[i][0], data[i][1], data[i][2], y[0], y[1], y[2]);
                            resy[j][i].back() = calc_derivative_backward(data[i][y.size() - 3], data[i][y.size() - 2], data[i].back(),
                                    y[y.size() - 3], y[y.size() - 2], y.back());
                    }
                    for (size_t i = 0; i < y.size(); ++i) {
                            resx[j][0][i] = calc_derivative_forward(data[0][i], data[1][i], data[2][i], x[0], x[1], x[2]);
                            resx[j].back()[i] = calc_derivative_backward(data[x.size() - 3][i], data[x.size() - 2][i], data.back()[i],
                                    x[x.size() - 3], x[x.size() - 2], x.back());
                    }
                    for (size_t i = 1; i < y.size() - 1; ++i) {
                            resy[j][0][i] = (data[0][i + 1] - data[0][i - 1]) / (y[i + 1] - y[i - 1]);
                            resy[j].back()[i] = (data.back()[i + 1] - data.back()[i - 1]) / (y[i + 1] - y[i - 1]);
                    }
                }

                return std::make_tuple(
                    utils::linear_interpolated_data_2d<T, V>(x, y, resx),
                    utils::linear_interpolated_data_2d<T, V>(x, y, resy)
                );
            }

            template<typename Arg, typename FX, typename FY, typename FT, typename FZ>
            auto compute(
                    const Arg& x0, const Arg& y0,
                    const Arg& l1, const size_t& nl,
                    const Arg& a0, const Arg& a1, const size_t& na,
                    const FX& fx, const FY& fy, const FT& ft, const FZ& fz, size_t& j, const size_t& nj) {
                const auto fs = dork::rk4<Arg>(Arg(0), l1, nl)(fx, fy, ft, fz);

                types::vector3d_t<Arg> rx(nj, types::vector2d_t<Arg>(na)),
                                       ry(nj, types::vector2d_t<Arg>(na));

                const auto mesh_l = utils::mesh_1d(Arg(0), l1, nl);
                const auto mesh_a = utils::mesh_1d(a0, a1, na);

                for (j = 0; j < nj; ++j)
                    for (size_t i = 0; i < na; ++i) {
                        rx[j][i].reserve(nl);
                        ry[j][i].reserve(nl);
                        auto solver = fs(x0, y0, std::cos(mesh_a[i]), std::sin(mesh_a[i]));
                        for (const auto& [l, x, y, t, z] : solver) {
                            rx[j][i].emplace_back(x);
                            ry[j][i].emplace_back(y);
                        }
                    }

                return 
                    std::make_tuple(
                        utils::linear_interpolated_data_2d<Arg>(mesh_a, mesh_l, std::move(rx)),
                        utils::linear_interpolated_data_2d<Arg>(mesh_a, mesh_l, std::move(ry)));
            }

        }// namespace __impl

        template<typename Arg>
        auto compute(
                const Arg& x0, const Arg& y0,
                const Arg& l1, const size_t& nl,
                const Arg& a0, const Arg& a1, const size_t& na,
                const utils::linear_interpolated_data_1d<Arg, Arg>& k_j) {
            types::vector1d_t<Arg> k0(k_j.size());
            for (size_t j = 0; j < k_j.size(); ++j)
                k0[j] = k_j[j].point(y0);

            const auto kd_j = __impl::calc_derivative(k_j);

            size_t j = 0;

            const auto fx = [&k0, &k_j, &j](const auto& s, const auto& x, const auto& y, const auto& t, const auto& z) {
                    return k0[j] * t / k_j[j].point(y);
            };

            const auto fy = [&k0, &k_j, &j](const auto& s, const auto& x, const auto& y, const auto& t, const auto& z) {
                    return k0[j] * z / k_j[j].point(y);
            };

            const auto ft = [&kd_j, &j](const auto& s, const auto& x, const auto& y, const auto& t, const auto& z) {
                    return 0;
            };

            const auto fz = [&kd_j, &j](const auto& s, const auto& x, const auto& y, const auto& t, const auto& z) {
                    return kd_j[j].point(y);
            };

            return __impl::compute(x0, y0, l1, nl, a0, a1, na, fx, fy, ft, fz, j, k_j.size());
        }

        template<typename Arg>
        auto compute(
                const Arg& x0, const Arg& y0,
                const Arg& l1, const size_t& nl,
                const Arg& a0, const Arg& a1, const size_t& na,
                const utils::linear_interpolated_data_2d<Arg, Arg>& k_j) {
            types::vector1d_t<Arg> k0(k_j.size());
            for (size_t j = 0; j < k_j.size(); ++j)
                k0[j] = k_j[j].point(x0, y0);

            const auto [kdx_j, kdy_j] = __impl::calc_derivatives(k_j);

            size_t j = 0;

            const auto fx = [&k0, &k_j, &j](const auto& s, const auto& x, const auto& y, const auto& t, const auto& z) {
                    return k0[j] * t / k_j[j].point(x, y);
            };

            const auto fy = [&k0, &k_j, &j](const auto& s, const auto& x, const auto& y, const auto& t, const auto& z) {
                    return k0[j] * z / k_j[j].point(x, y);
            };

            const auto ft = [&kdx_j, &j](const auto& s, const auto& x, const auto& y, const auto& t, const auto& z) {
                    return kdx_j[j].point(x, y);
            };

            const auto fz = [&kdy_j, &j](const auto& s, const auto& x, const auto& y, const auto& t, const auto& z) {
                    return kdy_j[j].point(x, y);
            };

            return __impl::compute(x0, y0, l1, nl, a0, a1, na, fx, fy, ft, fz, j, k_j.size());
        }

    }// namespace rays

}// namespace acstc