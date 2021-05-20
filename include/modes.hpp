#pragma once
#include <tuple>
#include <mutex>
#include <thread>
#include <cstddef>
#include <istream>
#include <algorithm>
#include <type_traits>
#include "normal_modes.h"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/assert.hpp"
#include "utils/callback.hpp"
#include "utils/interpolation.hpp"

namespace ample {

    template<typename T>
    class config;

    namespace _impl {

        template<typename T, typename V>
        struct modes_copier {

            static void copy(const NormalModes& n_m, types::vector3d_t<V>& k_j, types::vector4d_t<T>& phi_j,
                    const size_t& i, const size_t& j) {
                for (size_t k = 0; k < std::min(n_m.khs.size(), k_j.size()); ++k) {
                    k_j[k][i][j] = V(n_m.khs[k], n_m.mattenuation[k]);
                    phi_j[k][i][j] = std::move(n_m.mfunctions_zr[k]);
                }
            }

            static void copy(const NormalModes& n_m, types::vector2d_t<V>& k_j, types::vector3d_t<T>& phi_j,
                             const size_t& i) {
                for (size_t k = 0; k < std::min(n_m.khs.size(), k_j.size()); ++k) {
                    k_j[k][i] = V(n_m.khs[k], n_m.mattenuation[k]);
                    phi_j[k][i] = std::move(n_m.mfunctions_zr[k]);
                }
            }

        };

        template<typename T>
        struct modes_copier<T, T> {

            static void copy(const NormalModes& n_m, types::vector3d_t<T>& k_j, types::vector4d_t<T>& phi_j,
                    const size_t& i, const size_t& j) {
                for (size_t k = 0; k < std::min(n_m.khs.size(), k_j.size()); ++k) {
                    k_j[k][i][j] = n_m.khs[k];
                    phi_j[k][i][j] = std::move(n_m.mfunctions_zr[k]);
                }
            }

            static void copy(const NormalModes& n_m, types::vector2d_t<T>& k_j, types::vector3d_t<T>& phi_j,
                             const size_t& i) {
                for (size_t k = 0; k < std::min(n_m.khs.size(), k_j.size()); ++k) {
                    k_j[k][i] = n_m.khs[k];
                    phi_j[k][i] = std::move(n_m.mfunctions_zr[k]);
                }
            }

        };

    }// namespace _impl

    template<typename T = types::real_t, typename V = T>
    class modes {

    private:

        static constexpr auto Complex = !std::is_same_v<T, V>;

    public:

        explicit modes(const config<T>& config) : _config(config) {
            _n_m.iModesSubset = _config.mode_subset();
            _n_m.ppm = static_cast<unsigned int>(_config.ppm());
            _n_m.ordRich = static_cast<unsigned int>(_config.ord_rich());
            _n_m.f = _config.f();
            _n_m.M_betas = _config.betas();
            _n_m.eigen_type = "alglib";

            _n_m.M_depths.resize(_config.n_layers() + _config.bottom_layers().size());
            if (!config.additive_depth())
                std::copy(_config.bottom_layers().begin(), _config.bottom_layers().end(), _n_m.M_depths.begin() + _config.n_layers());

            _n_m.M_c1s.resize(_config.n_layers());
            _n_m.M_c1s.insert(_n_m.M_c1s.end(), _config.bottom_c1s().begin(), _config.bottom_c1s().end());

            _n_m.M_c2s.resize(_config.n_layers());
            _n_m.M_c2s.insert(_n_m.M_c2s.end(), _config.bottom_c2s().begin(), _config.bottom_c2s().end());

            _n_m.M_rhos.resize(config.n_layers(), T(1));
            _n_m.M_rhos.insert(_n_m.M_rhos.end(), _config.bottom_rhos().begin(), _config.bottom_rhos().end());

            _n_m.M_Ns_points.resize(_n_m.M_depths.size());
        }

        modes(const config<T>& config, const types::vector1d_t<T>& z) : modes(config) {
            _check_z(z);
            _n_m.zr.assign(z.begin(), z.end());
        }

        auto point(const T& x, const T& y, const size_t& c = -1) {
            _point(x, y, c);

            if constexpr (Complex) {
                types::vector1d_t<V> k_j(_n_m.khs.size());
                for (size_t j = 0; j < k_j.size(); ++j)
                    k_j[j] = V(_n_m.khs[j], _n_m.mattenuation[j]);
                return std::make_tuple(std::move(k_j), std::move(_n_m.mfunctions_zr));
            } else
                return std::make_tuple(std::move(_n_m.khs), std::move(_n_m.mfunctions_zr));            
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&>>>
        void line(const T& x, const T& y0, const T& y1, const size_t& ny, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto hy = (y1 - y0) / (ny - 1);
            const auto depth = _config.bathymetry().line(x, y0, y1, ny);

            _compute(ny, num_workers, [&](const size_t i0, const size_t i1, NormalModes n_m) {
                    auto y = y0 + hy * i0;
                    for (size_t i = i0; i < i1; ++i, y += hy) {
                        _point(n_m, x, y, depth[i], c);
                        callback(std::as_const(n_m), std::as_const(i));
                    }
                }
            );

        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&>>>
        void line(const T& x, const size_t& ny, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto [y0, y1] = _config.y_bounds();
            return line(x, y0, y1, ny, callback, num_workers, c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&>>>
        void line(const T& x, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            return line(x, _config.mny(), callback, num_workers, c);
        }

        auto vector_line(const T& x, const T& y0, const T& y1, const size_t& ny, const size_t& num_workers = 1, const size_t& c = -1) {
            return vector_line(x, y0, y1, ny, utils::nothing_callback(), num_workers, c);
        }

        auto vector_line(const T& x, const size_t& ny, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto [y0, y1] = _config.y_bounds();
            return vector_line(x, y0, y1, ny, num_workers, c);
        }

        auto vector_line(const T& x, const size_t& c = -1) {
            return vector_line(x, _config.mny(), c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&>>>
        auto vector_line(const T& x, const T& y0, const T& y1, const size_t& ny, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            size_t m = 0;
            types::vector2d_t<V> k_j;
            types::vector3d_t<T> phi_j;

            std::mutex mutex;
            line(x, y0, y1, ny, 
                utils::callbacks(
                    callback,
                    [&, mm=_config.max_mode()](const NormalModes& n_m, const size_t& i) mutable {
                        std::lock_guard<std::mutex> lock(mutex);
                        _fill_data(n_m, k_j, phi_j, ny, i, mm, m);
                    }
                ), num_workers, c
            );

            return std::make_tuple(std::move(k_j), std::move(phi_j));
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&>>>
        auto vector_line(const T& x, const size_t& ny, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto [y0, y1] = _config.y_bounds();
            return vector_line(x, y0, y1, ny, callback, num_workers, c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&>>>
        auto vector_line(const T& x, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            return vector_line(x, callback, num_workers, c);
        }

        auto interpolated_line(const T& x, const T& y0, const T& y1, const size_t& ny, const size_t& num_workers = 1, const size_t& c = -1) {
            return interpolated_line(x, y0, y1, ny, utils::nothing_callback(), num_workers, c);
        }

        auto interpolated_line(const T& x, const size_t& ny, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto [y0, y1] = _config.y_bounds();
            return interpolated_line(x, y0, y1, ny, num_workers, c);
        }

        auto interpolated_line(const T& x, const size_t& num_workers = 1, const size_t& c = -1) {
            return interpolated_line(x, _config.mny(), num_workers, c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&>>>
        auto interpolated_line(const T& x, const T& y0, const T& y1, const size_t& ny, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto y = utils::mesh_1d(y0, y1, ny);
            auto [k_j, phi_j] = vector_line(x, y0, y1, ny, callback, num_workers, c);

            return std::make_tuple(
                utils::linear_interpolated_data_1d<T, V>(y, std::move(k_j)),
                utils::linear_interpolated_data_2d<T, T>(y, _n_m.zr, std::move(phi_j))
            );
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&>>>
        auto interpolated_line(const T& x, const size_t& ny, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto [y0, y1] = _config.y_bounds();
            return interpolated_line(x, y0, y1, ny, callback, num_workers, c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&>>>
        auto interpolated_line(const T& x, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            return interpolated_line(x, _config.mny(), callback, num_workers, c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&, const size_t&>>>
        void field(
            const T& x0, const T& x1, const size_t& nx,
            const T& y0, const T& y1, const size_t& ny,
            C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto hx = (x1 - x0) / (nx - 1);
            const auto hy = (y1 - y0) / (ny - 1);
            const auto depth = _config.bathymetry().field(x0, x1, nx, y0, y1, ny);

            _compute(ny, num_workers, [&](const size_t j0, const size_t j1, NormalModes n_m) {
                    auto x = x0;
                    for (size_t i = 0; i < nx; ++i, x += hx) {
                        auto y = y0 + j0 * hy;
                        for (size_t j = j0; j < j1; ++j, y += hy) {
                            _point(n_m, x, y, depth[i][j], c);
                            callback(std::as_const(n_m), std::as_const(i), std::as_const(j));
                        }
                    }
                }
            );
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&, const size_t&>>>
        void field(const size_t& nx, const size_t& ny, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto [x0, x1] = _config.x_bounds();
            const auto [y0, y1] = _config.y_bounds();
            return field(x0, x1, nx, y0, y1, ny, callback, num_workers, c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&, const size_t&>>>
        void field(C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            return field(_config.mnx(), _config.mny(), callback, num_workers, c);
        }

        auto vector_field(
            const T& x0, const T& x1, const size_t& nx,
            const T& y0, const T& y1, const size_t& ny,
            const size_t& num_workers = 1, const size_t& c = -1) {
            return vector_field(x0, x1, nx, y0, y1, ny, utils::nothing_callback(), num_workers, c);
        }

        auto vector_field(const size_t& nx, const size_t& ny, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto [x0, x1] = _config.x_bounds();
            const auto [y0, y1] = _config.y_bounds();
            return vector_field(x0, x1, nx, y0, y1, ny, utils::nothing_callback(), num_workers, c);
        }

        auto vector_field(const size_t& num_workers = 1, const size_t& c = -1) {
            return vector_field(_config.mnx(), _config.mny(), utils::nothing_callback(), num_workers, c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&, const size_t&>>>
        auto vector_field(
            const T& x0, const T& x1, const size_t& nx,
            const T& y0, const T& y1, const size_t& ny,
            C&& callback, const size_t& num_workers = 1, const size_t& c = -1, const bool& smooth = false) {
            types::vector3d_t<V> k_j;
            types::vector4d_t<T> phi_j;

            std::mutex mutex;

            field(x0, x1, nx, y0, y1, ny,
                utils::callbacks(
                    callback,
                    [&, m=size_t(0), mm=_config.max_mode()](const NormalModes& n_m, const size_t& i, const size_t& j) mutable {
                        std::lock_guard<std::mutex> lock(mutex);
                        m = j ? m : 0;
                        _fill_data(n_m, k_j, phi_j, nx, ny, i, j, mm, m);
                    }
                ), num_workers, c
            );

            if (smooth)
                _smooth((y1 - y0) / (ny - 1), _config.border_width(), k_j, phi_j);

            return std::make_tuple(std::move(k_j), std::move(phi_j));
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&, const size_t&>>>
        auto vector_field(const size_t& nx, const size_t& ny, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto [x0, x1] = _config.x_bounds();
            const auto [y0, y1] = _config.y_bounds();
            return vector_field(x0, x1, nx, y0, y1, ny, callback, num_workers, c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&, const size_t&>>>
        auto vector_field(C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            return vector_field(_config.mnx(), _config.mny(), callback, num_workers, c);
        }

        auto interpolated_field(
            const T& x0, const T& x1, const size_t& nx,
            const T& y0, const T& y1, const size_t& ny,
            const size_t& num_workers = 1, const size_t& c = -1) {
            return interpolated_field(x0, x1, nx, y0, y1, ny, utils::nothing_callback(), num_workers, c);
        }

        auto interpolated_field(const size_t& nx, const size_t& ny, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto [x0, x1] = _config.x_bounds();
            const auto [y0, y1] = _config.y_bounds();
            return interpolated_field(x0, x1, nx, y0, y1, ny, num_workers, c);
        }

        auto interpolated_field(const size_t& num_workers = 1, const size_t& c = -1) {
            return interpolated_field(_config.mnx(), _config.mny(), num_workers, c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&, const size_t&>>>
        auto interpolated_field(
            const T& x0, const T& x1, const size_t& nx,
            const T& y0, const T& y1, const size_t& ny,
            C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto x = utils::mesh_1d(x0, x1, nx);
            const auto y = utils::mesh_1d(y0, y1, ny);
            auto [k_j, phi_j] = vector_field(x0, x1, nx, y0, y1, ny, callback, num_workers, c);

            return std::make_tuple(
                utils::linear_interpolated_data_2d<T, V>(x, y, std::move(k_j)),
                utils::linear_interpolated_data_3d<T, T>(x, y, _n_m.zr, std::move(phi_j))
            );
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&, const size_t&>>>
        auto interpolated_field(const size_t& nx, const size_t& ny, C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            const auto [x0, x1] = _config.x_bounds();
            const auto [y0, y1] = _config.y_bounds();
            return interpolated_field(x0, x1, nx, y0, y1, ny, callback, num_workers, c);
        }

        template<typename C, typename = std::enable_if_t<std::is_invocable_v<C, const NormalModes&, const size_t&, const size_t&>>>
        auto interpolated_field(C&& callback, const size_t& num_workers = 1, const size_t& c = -1) {
            return interpolated_field(_config.mnx(), _config.mny(), callback, num_workers, c);
        }

        void set_z(const types::vector1d_t<T>& z) {
            _check_z(z);
            _n_m.zr.assign(z.begin(), z.end());
        }

    private:

        NormalModes _n_m;
        const config<T>& _config;

        static constexpr T eps = 1;

        template<typename C>
        void _compute(const size_t& n, size_t num_workers, C&& callback) {
            num_workers = std::min(num_workers, n);

            if (num_workers <= 1) {
                callback(0, n, _n_m);
                return;
            }

            types::vector1d_t<std::thread> workers;
            workers.reserve(num_workers);

            const auto m = n / num_workers;

            for (size_t i = 0; i < num_workers; ++i)
                workers.emplace_back(
                        [&callback](const size_t& i, const size_t& j, NormalModes n_m) { callback(i, j, std::move(n_m)); },
                        m * i, i == num_workers - 1 ? n : m * (i + 1), _n_m);

            for (auto& it : workers)
                it.join();
        }

        static void _check_z(const types::vector1d_t<T>& z) {
            utils::dynamic_assert(std::all_of(z.begin(), z.end(), [](const auto& v) { return v >= -eps; }),
                  "All z coordinates must be non-negative");
        }

        static void _smooth(const T& h, const size_t& count, types::vector3d_t<V>& k_j, types::vector4d_t<T>& phi_j) {
            types::vector1d_t<T> coefficients(count);
            const auto d = (count - 1) * h;
            for (size_t i = 0; i < count; ++i)
                coefficients[i] = i * h / d;

            for (size_t j = 0; j < k_j.size(); ++j) {
                const auto& k0 = k_j[j][0];
                const auto& phi0 = phi_j[j][0];
                for (size_t i = 1; i < k_j[j].size(); ++i)
                    for (size_t l = 0, r = k0.size() - count; l < count; ++l, ++r) {
                        k_j[j][i][l] = k0[l] * (T(1) - coefficients[l]) + coefficients[l] * k_j[j][i][l];
                        k_j[j][i][r] = k_j[j][i][r] * (T(1) - coefficients[l]) + coefficients[l] * k0[r];
                        for (size_t k = 0; k < phi0[l].size(); ++k) {
                            phi_j[j][i][l][k] = phi0[l][k] * (T(1) - coefficients[l]) + coefficients[l] * phi_j[j][i][l][k];
                            phi_j[j][i][r][k] = phi_j[j][i][r][k] * (T(1) - coefficients[l]) + coefficients[l] * phi0[r][k];
                        }
                    }
            }
        }

        void _point(NormalModes& n_m, const T& x, const T& y, const T& depth, const size_t& c = -1) {
            utils::dynamic_assert(!n_m.zr.empty(), "There must be at least one depth value");

            if (depth <= eps) {
                n_m.khs.clear();
                n_m.mfunctions_zr.clear();
                return;
            }

            n_m.nmod = static_cast<int>(c == -1 ? _config.n_modes() : c);
            n_m.alpha = M_PI / 180 * (n_m.nmod > 0);

            auto buff = utils::mesh_1d(T(0), depth, _config.n_layers() + 1);
            std::copy(buff.begin() + 1, buff.end(), n_m.M_depths.begin());

            if (_config.additive_depth())
                std::transform(_config.bottom_layers().begin(), _config.bottom_layers().end(),
                               n_m.M_depths.begin() + _config.n_layers(), [&depth](const auto& z) { return z + depth; });

            _config.hydrology().line(x, T(0), depth, buff);
            std::copy(buff.begin(), buff.end() - 1, n_m.M_c1s.begin());
            std::copy(buff.begin() + 1, buff.end(), n_m.M_c2s.begin());

            n_m.M_Ns_points[0] = static_cast<unsigned>(std::round(n_m.ppm * n_m.M_depths[0]));
            for (size_t i = 1; i < n_m.M_depths.size(); ++i)
                n_m.M_Ns_points[i] = static_cast<unsigned>(std::round(n_m.ppm * (n_m.M_depths[i] - n_m.M_depths[i - 1])));

            n_m.compute_khs();
            n_m.compute_mfunctions_zr();
            if constexpr (Complex)
                n_m.compute_mattenuation();
        }

        void _point(const T& x, const T& y, const T& depth, const size_t& c = -1) {
            _point(_n_m, x, y, depth, c);
        }

        void _point(const T& x, const T& y, const size_t& c = -1) {
            _point(x, y, _config.bathymetry().point(x, y), c);
        }

        static auto _fill_data(const NormalModes& n_m, types::vector3d_t<V>& k_j, types::vector4d_t<T>& phi_j,
                const size_t& nx, const size_t& ny, const size_t& i, const size_t& j, const size_t& mm, size_t& m) {
            const auto n = std::min(n_m.khs.size(), mm);
            if (n > k_j.size()) {
                k_j.resize(n, types::vector2d_t<V>(nx, types::vector1d_t<V>(ny, V(0))));
                phi_j.resize(n, types::vector3d_t<T>(nx, types::vector2d_t<T>(ny, types::vector1d_t<T>(n_m.zr.size(), T(0)))));
            }
            _impl::modes_copier<T, V>::copy(n_m, k_j, phi_j, i, j);
            for (size_t k = m; k < n; ++k)
                for (size_t l = 0; l < j; ++l)
                    k_j[k][i][l] = k_j[k][i][j];
            for (size_t k = n; k < m; ++k)
                k_j[k][i][j] = k_j[k][i][j - 1];
            m = std::max(n, m);
        }

        static auto _fill_data(const NormalModes& n_m, types::vector2d_t<V>& k_j, types::vector3d_t<T>& phi_j,
                               const size_t& ny, const size_t& i, const size_t& mm, size_t& m) {
            const auto n = std::min(n_m.khs.size(), mm);
            if (n > k_j.size()) {
                k_j.resize(n, types::vector1d_t<V>(ny, V(0)));
                phi_j.resize(n, types::vector2d_t<T>(ny, types::vector1d_t<T>(n_m.zr.size(), T(0))));
            }
            _impl::modes_copier<T, V>::copy(n_m, k_j, phi_j, i);
            for (size_t k = m; k < n; ++k)
                for (size_t l = 0; l < i; ++l)
                    k_j[k][l] = k_j[k][i];
            for (size_t k = n; k < m; ++k)
                k_j[k][i] = k_j[k][i - 1];
            m = std::max(n, m);
        }

    };

}// namespace ample
