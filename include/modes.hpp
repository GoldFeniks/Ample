#pragma once
#include <tuple>
#include <cstddef>
#include <istream>
#include <algorithm>
#include "normal_modes.h"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/progress_bar.hpp"
#include "utils/interpolation.hpp"

namespace acstc {

    template<typename T>
    class config;

    template<typename T, typename V>
    class modes;

    namespace __impl {

        template<typename T>
        auto read_coords(std::istream& stream, types::vector1d_t<T>& v) {
            for (size_t i = 0; i < v.size(); ++i)
                stream >> v[i];
        }

        template<typename T, typename V>
        struct modes_reader {

            static auto read(std::istream& stream, types::vector3d_t<V>& res) {
                for (size_t i = 0; i < res.size(); ++i)
                    read(stream, res[i]);
            }

            static auto read(std::istream& stream, types::vector2d_t<V>& res) {
                T a, b;
                for (size_t i = 0; i < res.size(); ++i)
                    for (size_t j = 0; j < res[i].size(); ++j) {
                        stream >> a >> b;
                        res[i][j] = V(a, b);
                    }
            }

        };

        template<typename T>
        struct modes_reader<T, T> {

            static auto read(std::istream& stream, types::vector3d_t<T>& res) {
                for (size_t i = 0; i < res.size(); ++i)
                    read(stream, res[i]);
            }

            static auto read(std::istream& stream, types::vector2d_t<T>& res) {
                for (size_t i = 0; i < res.size(); ++i)
                    for (size_t j = 0; j < res[i].size(); ++j)
                        stream >> res[i][j];
            }

        };

        template<typename T, typename V>
        auto from_text(std::istream& stream, const size_t count) {
            size_t n, m, k;
            stream >> n >> m >> k;
            types::vector1d_t<T> x(n), y(m);
            types::vector3d_t<V> k_j(k, types::vector2d_t<V>(n, types::vector1d_t<V>(m)));
            types::vector3d_t<T> phi_j(k, types::vector2d_t<T>(n, types::vector1d_t<T>(m)));
            read_coords(stream, x);
            read_coords(stream, y);
            modes_reader<T, V>::read(stream, k_j);
            modes_reader<T, T>::read(stream, phi_j);
            ::acstc::modes<T, V>::smooth((y.back() - y.front()) / (y.size() - 1), count, k_j, phi_j);
            return std::make_tuple(
                    utils::linear_interpolated_data_2d<T, V>(x, y, std::move(k_j)),
                    utils::linear_interpolated_data_2d<T, T>(x, y, std::move(phi_j)));
        }

        template<typename T, typename V>
        auto const_from_text(std::istream& stream) {
            size_t n, k;
            stream >> n >> k;
            types::vector1d_t<T> y(n);
            types::vector2d_t<V> k_j(k, types::vector1d_t<V>(n));
            types::vector2d_t<T> phi_j(k, types::vector1d_t<T>(n));
            read_coords(stream, y);
            modes_reader<T, V>::read(stream, k_j);
            modes_reader<T, T>::read(stream, phi_j);
            return std::make_tuple(
                    utils::linear_interpolated_data_1d<T, V>(y, std::move(k_j)),
                    utils::linear_interpolated_data_1d<T, T>(y, std::move(phi_j)));
        }

        template<typename T, typename V>
        struct modes_copier {

            static void copy(const NormalModes& n_m, types::vector3d_t<V>& k_j, types::vector3d_t<T>& phi_j,
                    const size_t& i, const size_t& j) {
                for (size_t k = 0; k < std::min(n_m.mattenuation.size(), k_j.size()); ++k) {
                    k_j[k][i][j] = V(n_m.khs[k], n_m.mattenuation[k]);
                    phi_j[k][i][j] = n_m.mfunctions_zr[k][0];
                }
            }

            static void copy(const NormalModes& n_m, types::vector2d_t<V>& k_j, types::vector2d_t<T>& phi_j,
                             const size_t& i) {
                for (size_t k = 0; k < std::min(n_m.mattenuation.size(), k_j.size()); ++k) {
                    k_j[k][i] = V(n_m.khs[k], n_m.mattenuation[k]);
                    phi_j[k][i] = n_m.mfunctions_zr[k][0];
                }
            }

        };

        template<typename T>
        struct modes_copier<T, T> {

            static void copy(const NormalModes& n_m, types::vector3d_t<T>& k_j, types::vector3d_t<T>& phi_j,
                             const size_t& i, const size_t& j) {
                for (size_t k = 0; k < std::min(n_m.khs.size(), k_j.size()); ++k) {
                    k_j[k][i][j] = n_m.khs[k];
                    phi_j[k][i][j] = n_m.mfunctions_zr[k][0];
                }
            }

            static void copy(const NormalModes& n_m, types::vector2d_t<T>& k_j, types::vector2d_t<T>& phi_j,
                             const size_t& i) {
                for (size_t k = 0; k < std::min(n_m.khs.size(), k_j.size()); ++k) {
                    k_j[k][i] = n_m.khs[k];
                    phi_j[k][i] = n_m.mfunctions_zr[k][0];
                }
            }

        };

    }// namespace __impl

    template<typename T = types::real_t, typename V = T>
    class modes {

    public:

        modes() = delete;

        static auto create(const config<T>& config, const T& z, const bool show_progress = false) {
            return _create(config, config.bathymetry().x(), config.bathymetry().y(), config.bathymetry().data(), z, show_progress);
        }

        static auto create(const config<T>& config, const T& z,
                const T& x0, const T& x1, const size_t& nx,
                const T& y0, const T& y1, const size_t& ny,
                const bool show_progress = false) {
            return _create(config, utils::mesh_1d(x0, x1, nx), utils::mesh_1d(y0, y1, ny),
                    config.bathymetry().field(x0, x1, nx, y0, y1, ny), z, show_progress);
        }

        static auto create(const config<T>& config, const T& z, const size_t& nx, const size_t& ny, const bool show_progress = false) {
            const auto [x0, x1, y0, y1] = config.bounds();
            return create(config, z, x0, x1, nx, y0, y1, ny, show_progress);
        }

        static auto create(const config<T>& config, const T& z, const T& y0, const T& y1, const size_t& ny, const bool show_progress = false) {
            return _create(config, utils::mesh_1d(y0, y1, ny),
                    config.bathymetry().line(config.bathymetry().x().front(), y0, y1, ny), z, show_progress);
        }

        static auto create(const config<T>& config, const T& z, const size_t& ny, const bool show_progress = false) {
            const auto [y0, y1] = config.y_bounds();
            return create(config, z, y0, y1, ny, show_progress);
        }

        static auto from_text(std::istream& stream, const size_t count) {
            return __impl::from_text<T, V>(stream, count);
        }

        static auto from_text(std::istream&& stream, const size_t count) {
            return from_text(stream, count);
        }

        template<typename S = uint32_t>
        static auto from_binary(std::istream& stream, const size_t count) {
            S n, m, k;
            stream.read(reinterpret_cast<char*>(&n), sizeof(S));
            stream.read(reinterpret_cast<char*>(&m), sizeof(S));
            stream.read(reinterpret_cast<char*>(&k), sizeof(S));
            types::vector1d_t<T> x(n), y(m);
            stream.read(reinterpret_cast<char*>(x.data()), sizeof(T) * n);
            stream.read(reinterpret_cast<char*>(y.data()), sizeof(T) * m);
            types::vector3d_t<V> k_j(k, types::vector2d_t<V>(n, types::vector1d_t<V>(m)));
            types::vector3d_t<T> phi_j(k, types::vector2d_t<T>(n, types::vector1d_t<T>(m)));
            for (size_t j = 0; j < k; ++j)
                for (size_t i = 0; i < n; ++i)
                    stream.read(reinterpret_cast<char*>(k_j[j][i].data()), sizeof(V) * m);
            for (size_t j = 0; j < k; ++j)
                for (size_t i = 0; i < n; ++i)
                    stream.read(reinterpret_cast<char*>(phi_j[j][i].data()), sizeof(V) * m);
            smooth((y.back() - y.front()) / (y.size() - 1), count, k_j, phi_j);
            return std::make_tuple(
                    utils::linear_interpolated_data_2d<T, V>(x, y, std::move(k_j)),
                    utils::linear_interpolated_data_2d<T, T>(x, y, std::move(phi_j)));
        }

        template<typename S = uint32_t>
        static auto from_binary(std::istream&& stream, const size_t count) {
            return from_binary<S>(stream, count);
        }

        static auto const_from_text(std::istream& stream) {
            return __impl::const_from_text<T, V>(stream);
        }

        static auto const_from_text(std::istream&& stream) {
            return const_from_text(stream);
        }

        template<typename S = uint32_t>
        static auto const_from_binary(std::istream& stream) {
            S n, k;
            stream.read(reinterpret_cast<char*>(&n), sizeof(S));
            stream.read(reinterpret_cast<char*>(&k), sizeof(S));
            types::vector1d_t<T> y(n);
            types::vector2d_t<V> k_j(k, types::vector1d_t<V>(n));
            types::vector2d_t<T> phi_j(k, types::vector1d_t<T>(n));
            for (size_t j = 0; j < k; ++j)
                stream.read(reinterpret_cast<char*>(k_j[j].data()), sizeof(V) * n);
            for (size_t j = 0; j < k; ++j)
                stream.read(reinterpret_cast<char*>(phi_j[j].data()), sizeof(V) * n);
            return std::make_tuple(
                    utils::linear_interpolated_data_1d<T, V>(y, std::move(k_j)),
                    utils::linear_interpolated_data_1d<T, T>(y, std::move(phi_j)));
        }

        template<typename S = uint32_t>
        static auto const_from_binary(std::istream&& stream) {
            return const_from_binary(stream);
        }

        static auto calc_modes(const config<T>& config, const T& x, const T& depth, const T& z) {
            NormalModes n_m;
            n_m.iModesSubset = config.mode_subset();
            n_m.ppm = static_cast<unsigned int>(config.ppm());
            n_m.ordRich = static_cast<unsigned int>(config.ordRich());
            n_m.zr.push_back(z);
            n_m.f = config.f();
            n_m.M_betas = config.betas();

            auto buff = utils::mesh_1d(T(0), depth, config.n_layers() + 1);
            n_m.M_depths.assign(buff.begin() + 1, buff.end());
            if (config.additive_depth()) {
                n_m.M_depths.reserve(n_m.M_depths.size() + config.bottom_layers().size());
                for (const auto& it : config.bottom_layers())
                    n_m.M_depths.push_back(it + depth);
            } else
                n_m.M_depths.insert(n_m.M_depths.end(), config.bottom_layers().begin(), config.bottom_layers().end());

            config.hydrology().line(x, T(0), depth, buff);
            n_m.M_c1s.assign(buff.begin(), buff.end() - 1);
            n_m.M_c1s.insert(n_m.M_c1s.end(), config.bottom_c1s().begin(), config.bottom_c1s().end());
            n_m.M_c2s.assign(buff.begin() + 1, buff.end());
            n_m.M_c2s.insert(n_m.M_c2s.end(), config.bottom_c2s().begin(), config.bottom_c2s().end());

            n_m.M_rhos.resize(config.n_layers(), T(1));
            n_m.M_rhos.insert(n_m.M_rhos.end(), config.bottom_rhos().begin(), config.bottom_rhos().end());

            n_m.M_Ns_points.resize(n_m.M_depths.size());
            n_m.M_Ns_points[0] = static_cast<unsigned>(std::round(round(n_m.ppm * n_m.M_depths[0])));
            for (size_t i = 1; i < n_m.M_depths.size(); ++i)
                n_m.M_Ns_points[i] = static_cast<unsigned>(round(n_m.ppm * (n_m.M_depths[i] - n_m.M_depths[i - 1])));

            n_m.compute_khs();
            n_m.compute_mfunctions_zr();
            if (config.complex_modes())
                n_m.compute_mattenuation();
            return n_m;
        }

        static void smooth(const T& h, const size_t count, types::vector3d_t<V>& k, types::vector3d_t<T>& phi) {
            types::vector1d_t<T> coefficients(count);
            const auto d = (count - 1) * h;
            for (size_t i = 0; i < count; ++i)
                coefficients[i] = i * h / d;
            for (size_t j = 0; j < k.size(); ++j) {
                const auto& k0 = k[j][0];
                const auto& phi0 = phi[j][0];
                for (size_t i = 1; i < k.size(); ++i)
                    for (size_t l = 0, r = k0.size() - count; l < count; ++l, ++r) {
                        k[j][i][l] = k0[l] * (T(1) - coefficients[l]) + coefficients[l] * k[j][i][l];
                        k[j][i][r] = k[j][i][r] * (T(1) - coefficients[l]) + coefficients[l] * k0[r];
                        phi[j][i][l] = phi0[l] * (T(1) - coefficients[l]) + coefficients[l] * phi[j][i][l];
                        phi[j][i][r] = phi[j][i][r] * (T(1) - coefficients[l]) + coefficients[l] * phi0[r];
                    }
            }
        }

    private:

        template<typename XV, typename YV, typename DV>
        static auto _create(const config<T>& config, const XV& x, const YV& y, const DV& data, const T& z, const bool show_progress) {
            types::vector3d_t<V> k_j;
            types::vector3d_t<T> phi_j;
            const auto mm = config.max_mode();
            const auto nx = x.size();
            const auto ny = y.size();
            utils::progress_bar pbar(nx * ny, "Modes");

            for (size_t i = 0; i < nx; ++i) {
                size_t m = 0;
                for (size_t j = 0; j < ny; ++j) {
                    _fill_data(calc_modes(config, x[i], data[i][j], z), k_j, phi_j, nx, ny, i, j, mm, m);
                    if (show_progress)
                        pbar();
                }
            }
            smooth((y.back() - y.front()) / (y.size() - 1), config.border_width(), k_j, phi_j);
            return std::make_tuple(
                    utils::linear_interpolated_data_2d<T, V>(x, y, std::move(k_j)),
                    utils::linear_interpolated_data_2d<T, T>(x, y, std::move(phi_j)));
        }

        template<typename YV, typename DV>
        static auto _create(const config<T>& config, const YV& y, const DV& data, const T& z, const bool show_progress) {
            types::vector2d_t<V> k_j;
            types::vector2d_t<T> phi_j;
            const auto mm = config.max_mode();
            const auto ny = y.size();
            utils::progress_bar pbar(ny, "Modes");

            size_t m = 0;
            for (size_t i = 0; i < ny; ++i) {
                _fill_data(calc_modes(config, config.x0(), data[i], z), k_j, phi_j, ny, i, mm, m);
                if (show_progress)
                    pbar();
            }
            return std::make_tuple(
                    utils::linear_interpolated_data_1d<T, V>(y, std::move(k_j)),
                    utils::linear_interpolated_data_1d<T, T>(y, std::move(phi_j)));
        }

        static auto _fill_data(const NormalModes& n_m, types::vector3d_t<V>& k_j, types::vector3d_t<T>& phi_j,
                const size_t& nx, const size_t& ny, const size_t& i, const size_t& j, const size_t& mm, size_t& m) {
            const auto n = std::min(n_m.khs.size(), mm);
            if (n > k_j.size()) {
                k_j.resize(n, types::vector2d_t<V>(nx, types::vector1d_t<V>(ny, V(0))));
                phi_j.resize(n, types::vector2d_t<T>(nx, types::vector1d_t<T>(ny, T(0))));
            }
            __impl::modes_copier<T, V>::copy(n_m, k_j, phi_j, i, j);
            for (size_t k = m; k < n; ++k)
                for (size_t l = 0; l < j; ++l)
                    k_j[k][i][l] = k_j[k][i][j];
            for (size_t k = n; k < m; ++k)
                k_j[k][i][j] = k_j[k][i][j - 1];
            m = std::max(n, m);
        }

        static auto _fill_data(const NormalModes& n_m, types::vector2d_t<V>& k_j, types::vector2d_t<T>& phi_j,
                               const size_t& ny, const size_t& i, const size_t& mm, size_t& m) {
            const auto n = std::min(n_m.khs.size(), mm);
            if (n > k_j.size()) {
                k_j.resize(n, types::vector1d_t<V>(ny, V(0)));
                phi_j.resize(n, types::vector1d_t<T>(ny, T(0)));
            }
            __impl::modes_copier<T, V>::copy(n_m, k_j, phi_j, i);
            for (size_t k = m; k < n; ++k)
                for (size_t l = 0; l < i; ++l)
                    k_j[k][l] = k_j[k][i];
            for (size_t k = n; k < m; ++k)
                k_j[k][i] = k_j[k][i - 1];
            m = std::max(n, m);
        }

    };

}// namespace acstc
