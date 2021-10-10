#pragma once

#include <cmath>
#include <cstddef>
#include <algorithm>
#include <type_traits>
#include "dork.hpp"
#include "utils/types.hpp"
#include "utils/progress_bar.hpp"
#include "utils/interpolation.hpp"

namespace ample::rays {

    namespace _impl {

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
        auto calc_derivative_x(const utils::linear_interpolated_data_2d<T, V>& values) {
            const auto& x = values.template get<0>();
            const auto& y = values.template get<1>();
            types::vector3d_t<V> resx(values.size());
            for (size_t j = 0; j < values.size(); ++j) {
                const auto& data = values[j].data();
                resx[j].resize(x.size(), types::vector1d_t<V>(y.size()));
                for (size_t i = 1; i < x.size() - 1; ++i) {
                    resx[j][i][0] = (data[i + 1][0] - data[i - 1][0]) / (x[i + 1] - x[i - 1]);
                    resx[j][i].back() = (data[i + 1].back() - data[i - 1].back()) / (x[i + 1] - x[i - 1]);
                    for (size_t k = 1; k < y.size() - 1; ++k)
                        resx[j][i][k] = (data[i + 1][k] - data[i - 1][k]) / (x[i + 1] - x[i - 1]);
                }
                for (size_t i = 0; i < y.size(); ++i) {
                    resx[j][0][i] = calc_derivative_forward(data[0][i], data[1][i], data[2][i], x[0], x[1], x[2]);
                    resx[j].back()[i] = calc_derivative_backward(data[x.size() - 3][i], data[x.size() - 2][i], data.back()[i],
                            x[x.size() - 3], x[x.size() - 2], x.back());
                }
            }
            return utils::linear_interpolated_data_2d<T, V>(x, y, resx);
        }

        template<typename T, typename V>
        auto calc_derivatives(const utils::linear_interpolated_data_2d<T, V>& values) {
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

        template<typename Arg, typename FX, typename FY, typename FT, typename FZ, typename C, typename = std::enable_if_t<std::is_invocable_v<C, size_t, size_t, size_t, Arg, Arg, Arg, Arg>>>
        void compute(
                const Arg& x0, const Arg& y0,
                const Arg& l1, const size_t& nl,
                const Arg& a0, const Arg& a1, const size_t& na,
                const FX& fx, const FY& fy, const FT& ft, const FZ& fz, size_t& j, const size_t& nj,
                C&& callback,
                const bool& show_progress) {
            const auto fs = dork::rk4<Arg>(Arg(0), l1, nl)(fx, fy, ft, fz);

            const auto mesh_a = utils::mesh_1d(a0, a1, na);

            utils::progress_bar pbar(nj * na * nl, "Rays", show_progress);

            size_t k = 0;
            for (j = 0; j < nj; ++j)
                for (size_t i = 0; i < na; ++i, ++k) {
                    auto solver = fs(x0, y0, std::cos(mesh_a[i]), std::sin(mesh_a[i]));
                    for (const auto& [l, x, y, t, z] : solver) {
                        callback(j, i, k, x, y, mesh_a[i], l);
                        pbar.next();
                    }
                }
        }

    }// namespace _impl

    template<typename Arg, typename Val, typename C, typename = std::enable_if_t<std::is_invocable_v<C, size_t, size_t, size_t, Arg, Arg, Arg, Arg>>>
    void compute(
            const Arg& x0, const Arg& y0,
            const Arg& l1, const size_t& nl,
            const Arg& a0, const Arg& a1, const size_t& na,
            const utils::linear_interpolated_data_1d<Arg, Val>& k_j,
            C&& callback,
            const bool& show_progress = false) {
        if constexpr (!std::is_same_v<Arg, Val>) {
            const auto& x = k_j.template get<0>();

            types::vector2d_t<Arg> data(k_j.size(), types::vector1d_t<Arg>(x.size()));
            for (size_t i = 0; i < data.size(); ++i) {
                const auto& k_j_data = k_j[i].data();
                std::transform(k_j_data.begin(), k_j_data.end(), data[i].begin(),
                    [](const auto& value) { return value.real(); }
                );
            }

            compute(x0, y0, l1, nl, a0, a1, na, utils::linear_interpolated_data_1d<Arg, Arg>(x, data), callback, show_progress);
        } else {
            types::vector1d_t<Arg> k0(k_j.size());
            for (size_t j = 0; j < k_j.size(); ++j)
                k0[j] = k_j[j].point(y0);

            const auto kd_j = _impl::calc_derivative(k_j);

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

            _impl::compute(x0, y0, l1, nl, a0, a1, na, fx, fy, ft, fz, j, k_j.size(), callback, show_progress);
        }
    }

    template<typename Arg, typename Val, typename C, typename = std::enable_if_t<std::is_invocable_v<C, size_t, size_t, size_t, Arg, Arg, Arg, Arg>>>
    auto compute(
            const Arg& x0, const Arg& y0,
            const Arg& l1, const size_t& nl,
            const Arg& a0, const Arg& a1, const size_t& na,
            const utils::linear_interpolated_data_2d<Arg, Val>& k_j,
            C&& callback,
            const bool& show_progress = false) {
        if constexpr (!std::is_same_v<Arg, Val>) {
            const auto& x = k_j.template get<0>();
            const auto& y = k_j.template get<1>();

            types::vector3d_t<Arg> data(k_j.size(), types::vector2d_t<Arg>(x.size(), types::vector1d_t<Arg>(y.size())));
            for (size_t i = 0; i < data.size(); ++i) {
                const auto& k_j_data = k_j[i].data();
                for (size_t j = 0; j < data[i].size(); ++j)
                    std::transform(k_j_data[j].begin(), k_j_data[j].end(), data[i][j].begin(),
                        [](const auto& value) { return value.real(); }
                    );
            }

            return compute(x0, y0, l1, nl, a0, a1, na, utils::linear_interpolated_data_2d<Arg, Arg>(x, y, data), callback, show_progress);
        } else {
            types::vector1d_t<Arg> k0(k_j.size());
            for (size_t j = 0; j < k_j.size(); ++j)
                k0[j] = k_j[j].point(x0, y0);

            const auto [kdx_j, kdy_j] = _impl::calc_derivatives(k_j);

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

            return _impl::compute(x0, y0, l1, nl, a0, a1, na, fx, fy, ft, fz, j, k_j.size(), callback, show_progress);
        }
    }

    template<typename Arg, typename KJ>
    auto compute(
            const Arg& x0, const Arg& y0,
            const Arg& l1, const size_t& nl,
            const Arg& a0, const Arg& a1, const size_t& na,
            const KJ& k_j,
            const bool& show_progress = false) {
        const auto nj = k_j.size();

        types::vector3d_t<Arg> rx(nj, types::vector2d_t<Arg>(na)),
                               ry(nj, types::vector2d_t<Arg>(na));

        for (size_t j = 0; j < nj; ++j)
            for (size_t i = 0; i < na; ++i) {
                rx[j][i].reserve(nl);
                ry[j][i].reserve(nl);
            }

        compute(x0, y0, l1, nl, a0, a1, na, k_j, [&](const auto& j, const auto& i, const auto& k, const auto& x, const auto& y, const auto& a, const auto& l) {
            rx[j][i].emplace_back(x);
            ry[j][i].emplace_back(y);
        }, show_progress);

        const auto mesh_a = utils::mesh_1d(a0, a1, na);
        const auto mesh_l = utils::mesh_1d(Arg(0), l1, nl);

        return std::make_tuple(
                utils::linear_interpolated_data_2d<Arg>(mesh_a, mesh_l, std::move(rx)),
                utils::linear_interpolated_data_2d<Arg>(mesh_a, mesh_l, std::move(ry)));
    }

}// namespace ample::rays
