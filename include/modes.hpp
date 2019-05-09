#pragma once
#include <tuple>
#include "config.hpp"
#include "normal_modes.h"
#include "utils/types.hpp"
#include "utils/interpolation.hpp"

namespace acstc {

    template<typename T = types::real_t, typename V = types::complex_t, typename I = utils::linear_interpolation>
    class modes {

    public:

        modes() = delete;

        template<typename CI>
        static auto create(const config<T, CI>& config) {
            return _create(config, config.bathymetry().x(), config.bathymetry().y(), config.bathymetry().data());
        }

        template<typename CI>
        static auto create(const config<T, CI>& config,
                const T& x0, const T& x1, const size_t& nx,
                const T& y0, const T& y1, const size_t& ny) {
            return _create(config, utils::mesh_1d(x0, x1, nx), utils::mesh_1d(y0, y1, ny),
                    config.bathymetry().field(x0, x1, nx, y0, y1, ny));
        }

        template<typename CI>
        static auto create(const config<T, CI>& config, const size_t& nx, const size_t& ny) {
            const auto [x0, x1, y0, y1] = config.bounds();
            return create(config, x0, x1, nx, y0, y1, ny);
        }

        template<typename CI>
        static auto create(const config<T, CI>& config, const T& y0, const T& y1, const size_t& ny) {
            return _create(config, utils::mesh_1d(y0, y1, ny),
                    config.bathymetry().line(config.bathymetry().x().front(), y0, y1, ny));
        }

        template<typename CI>
        static auto create(const config<T, CI>& config, const size_t& ny) {
            const auto [y0, y1] = config.y_bounds();
            return create(config, y0, y1, ny);
        }

    private:

        template<typename CI, typename XV, typename YV, typename DV>
        static auto _create(const config<T, CI>& config, const XV& x, const YV& y, const DV& data) {
            types::vector3d_t<V> k_j;
            types::vector3d_t<T> phi_j;
            const auto nx = x.size();
            const auto ny = y.size();
            for (size_t i = 0; i < nx; ++i) {
                size_t m = 0;
                for (size_t j = 0; j < ny; ++j)
                    _fill_data(_calc_modes(config, data[i][j]), k_j, phi_j, nx, ny, i, j, m);
            }
            return std::make_tuple(
                    utils::interpolated_data_2d(x, y, std::move(k_j)),
                    utils::interpolated_data_2d(x, y, std::move(phi_j)));
        }

        template<typename CI, typename YV, typename DV>
        static auto _create(const config<T, CI>& config, const YV& y, const DV& data) {
            types::vector2d_t<V> k_j;
            types::vector2d_t<T> phi_j;
            const auto ny = y.size();
            size_t m = 0;
            for (size_t i = 0; i < ny; ++i)
                _fill_data(_calc_modes(config, data[i]), k_j, phi_j, ny, i, m);
            return std::make_tuple(
                    utils::interpolated_data_1d(y, std::move(k_j)),
                    utils::interpolated_data_1d(y, std::move(phi_j)));
        }

        template<typename CI>
        static auto _calc_modes(const config<T, CI>& config, const T& depth) {
            NormalModes n_m;
            n_m.iModesSubset = config.mode_subset();
            n_m.ppm = static_cast<unsigned int>(config.ppm());
            n_m.ordRich = static_cast<unsigned int>(config.ordRich());
            n_m.zr.push_back(config.zr());
            n_m.f = config.f();
            n_m.M_c1s = { T(1500), T(1700) };
            n_m.M_c2s = { T(1500), T(1700) };
            n_m.M_rhos = { T(1), T(2) };
            n_m.M_depths = { depth, T(500) };
            n_m.M_Ns_points.resize(n_m.M_depths.size());
            n_m.M_Ns_points[0] = static_cast<unsigned>(std::round(round(n_m.ppm * n_m.M_depths[0])));
            for (size_t i = 1; i < n_m.M_depths.size(); ++i)
                n_m.M_Ns_points[i] = static_cast<unsigned>(round(n_m.ppm * (n_m.M_depths[i] - n_m.M_depths[i - 1])));
            n_m.compute_khs();
            n_m.compute_mfunctions_zr();
            return n_m;
        }

        static auto _fill_data(const NormalModes& n_m, types::vector3d_t<V>& k_j, types::vector3d_t<V>& phi_j,
                const size_t& nx, const size_t& ny, const size_t& i, const size_t& j, size_t& m) {
            const auto n = n_m.khs.size();
            if (n > k_j.size()) {
                k_j.resize(n, types::vector2d_t<V>(nx, types::vector1d_t<V>(ny, V(0))));
                phi_j.resize(n, types::vector2d_t<V>(nx, types::vector1d_t<V>(ny, V(0))));
            }
            for (size_t k = 0; k < n; ++k) {
                k_j[k][i][j] = n_m.khs[k];
                phi_j[k][i][j] = n_m.mfunctions_zr[0][k];
            }
            for (size_t l = 0; l < j; ++l)
                for (size_t k = m; k < n; ++k)
                    k_j[k][i][l] = k_j[k][i][j];
            for (size_t k = n; k < m; ++k)
                k_j[k][i][j] = k_j[k][i][j - 1];
            m = n;
        }

        static auto _fill_data(const NormalModes& n_m, types::vector2d_t<V>& k_j, types::vector2d_t<V>& phi_j,
                               const size_t& ny, const size_t& i, size_t& m) {
            const auto n = n_m.khs.size();
            if (n > k_j.size()) {
                k_j.resize(n, types::vector1d_t<V>(ny, V(0)));
                phi_j.resize(n, types::vector1d_t<V>(ny, V(0)));
            }
            for (size_t k = 0; k < n; ++k) {
                k_j[k][i] = n_m.khs[k];
                phi_j[k][i] = n_m.mfunctions_zr[0][k];
            }
            for (size_t l = 0; l < i; ++l)
                for (size_t k = m; k < n; ++k)
                    k_j[k][l] = k_j[k][i];
            for (size_t k = n; k < m; ++k)
                k_j[k][i] = k_j[k][i - 1];
            m = n;
        }

    };

}// namespace acstc
