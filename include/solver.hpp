#pragma once
#include <cmath>
#include <mutex>
#include <thread>
#include <cstddef>
#include <complex>
#include <algorithm>
#include <condition_variable>
#include "config.hpp"
#include "utils/event.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/assert.hpp"
#include "coefficients.hpp"
#include "threads/pool.hpp"
#include "boundary_conditions.hpp"
#include "utils/interpolation.hpp"
#include "threads/buffer_manager.hpp"

//WINDOWS DEFINES VERY NICE MACRO FOR NO REASON
#ifdef IN
#undef IN
#endif

namespace ample {
    
    using namespace std::complex_literals;

    template<typename BC, typename Arg = typename BC::arg_t, typename Val = typename BC::val_t>
    class solver {

    public:

        utils::event<const Arg&, const types::vector2d_t<Val>&> on_solution;
        utils::event<const size_t&, const Arg&, const Val&, const types::vector2d_t<Val>&> on_mode_solution;

        solver(BC boundary_conditions, coefficients<Val> coefficients,
               const Arg& x0, const Arg& x1, const size_t& nx,
               const Arg& y0, const Arg& y1, const size_t& ny,
               const Arg& z0, const Arg& z1, const size_t& nz)
                   : _nx(nx), _ny(ny), _nz(nz),
                     _boundary_conditions(std::move(boundary_conditions)),
                     _coefficients(std::move(coefficients)),
                     _hx((x1 - x0) / (nx - 1)),
                     _x0(x0), _y0(y0), _y1(y1), _z0(z0), _z1(z1) {}

        solver(BC boundary_conditions, const config<Arg>& config)
                  : solver(std::move(boundary_conditions), config.template make_coefficients<Val>(),
                           config.x0(), config.x1(), config.nx(),
                           config.y0(), config.y1(), config.ny(),
                           config.z0(), config.z1(), config.nz()) {}

        template<typename IN, typename VL>
        void solve(const IN& init,
                   const types::vector1d_t<VL>& k0,
                   const utils::linear_interpolated_data_1d<Arg, VL>& k_int,
                   const utils::linear_interpolated_data_2d<Arg, Arg>& phi_int,
                   const size_t num_workers = 1,
                   const size_t buff_size = 100) const {
            const auto nm = k0.size();

            if (nm == 0)
                return;

            utils::dynamic_assert(k_int.size() == nm && phi_int.size() == nm,
                "Inputs k0(", nm, "), k_int(", k_int.size(), "), phi_int(", phi_int.size(), ") must have the same size");

            types::vector1d_t<Val> a0(nm);
            types::vector2d_t<VL>  kk(nm);
            types::vector2d_t<Val> aa(nm), bb(nm);
            types::vector3d_t<Arg> ph(nm);
            for (size_t j = 0; j < nm; ++j) {
                kk[j] = k_int[j].line(_y0, _y1, _ny);
                ph[j] = phi_int[j].field(_y0, _y1, _ny, _z0, _z1, _nz);
                std::tie(a0[j], aa[j], bb[j]) = _coefficients.get(im * k0[j] * _hx);
            }

            const auto nw = _boundary_conditions.width();
            const auto nc = _coefficients.nc();

            auto band_builder = _boundary_conditions.get_band_builder(k0, bb, _y0, _y1, _ny);
            const auto [ac, bc, cc] = band_builder.coefficients();
            const auto ny = band_builder.ny();
            band_builder.update(kk);

            types::vector2d_t<Val> bv(_ny, types::vector1d_t<Val>(_nz, Val(0)));

            auto cv = init.make(band_builder.y0(), band_builder.y1(), ny, k0.size());

            for (size_t j = 0; j < nm; ++j) {
                const auto exp = std::exp(im * k0[j] * _x0);
                for (size_t y = 0, i = nw; y < _ny; ++y, ++i)
                    for (size_t z = 0; z < _nz; ++z)
                        bv[y][z] += ph[j][y][z] * cv[j][i] * exp;
            }

            on_solution(_x0, bv);

            auto solve_func = [&, &ac=ac, &bc=bc, &cc=cc](const size_t& j, const Arg& x, auto& ov, auto& nv, auto& solver) {
                nv.assign(nc, cv[j]);
                for (auto& it : nv)
                    it.front() = it.back() = ze;

                cv[j].front() = cv[j].back() = ze;
                for (size_t y = 1; y < ny - 1; ++y)
                    cv[j][y] *= a0[j];

                for (size_t i = 0; i < nc; ++i) {
                    solver(ac[j][i], bc[j][i], cc[j][i], nv[i]);
                    std::transform(cv[j].begin(), cv[j].end(), nv[i].begin(), cv[j].begin(),
                                   [&c=aa[j][i]](const auto& a, const auto& b) { return a + c * b; });
                }

                _complete_mode_solution(j, x, nw, k0[j], ov, ph[j], cv[j]);
            };

            _compute(solve_func, nm, ny, nc, num_workers, buff_size);
        }

        template<typename IN, typename VL>
        void solve(const IN& init,
                   const types::vector1d_t<VL>& k0,
                   const utils::linear_interpolated_data_2d<Arg, VL>& k_int,
                   const utils::linear_interpolated_data_3d<Arg, Arg>& phi_int,
                   const size_t num_workers = 1,
                   const size_t buff_size = 100) const {
            const auto nm = k0.size();

            if (nm == 0)
                return;

            utils::dynamic_assert(k_int.size() == nm && phi_int.size() == nm,
                "Inputs k0(", nm, "), k_int(", k_int.size(), "), phi_int(", phi_int.size(), ") must have the same size");

            types::vector1d_t<Val> a0(nm);
            types::vector2d_t<Val> aa(nm), bb(nm);
            for (size_t j = 0; j < nm; ++j)
                std::tie(a0[j], aa[j], bb[j]) = _coefficients.get(im * k0[j] * _hx);

            const auto nw = _boundary_conditions.width();
            const auto nc = _coefficients.nc();

            auto band_builder = _boundary_conditions.get_band_builder(k0, bb, _y0, _y1, _ny);
            const auto [ac, bc, cc] = band_builder.coefficients();
            const auto ny = band_builder.ny();

            types::vector2d_t<Arg> ip(_ny, types::vector1d_t<Arg>(_nz));
            types::vector2d_t<Val> bv(_ny, types::vector1d_t<Val>(_nz, Val(0)));

            auto cv = init.make(band_builder.y0(), band_builder.y1(), ny, k0.size());

            for (size_t j = 0; j < nm; ++j) {
                const auto exp = std::exp(im * k0[j] * _x0);

                phi_int[j].field(_x0, _y0, _y1, _z0, _z1, ip);
                for (size_t y = 0, i = nw; y < _ny; ++y, ++i)
                    for (size_t z = 0; z < _nz; ++z)
                        bv[y][z] += ip[y][z] * cv[j][i] * exp;
            }

            on_solution(_x0, bv);

            types::vector2d_t<VL>  kk(nm, types::vector1d_t<VL> (_ny));
            types::vector3d_t<Arg> ph(nm, types::vector2d_t<Arg>(_ny, types::vector1d_t<Arg>(_nz)));

            auto solve_func = [&, &ac=ac, &bc=bc, &cc=cc](const size_t& j, const Arg& x, auto& ov, auto& nv, auto& solver) {
                k_int[j].line(x, _y0, _y1, kk[j]);
                phi_int[j].field(x, _y0, _y1, _z0, _z1, ph[j]);

                band_builder.update(kk, j, j + 1);

                nv.assign(nc, cv[j]);
                for (auto &it : nv)
                    it.front() = it.back() = ze;

                cv[j].front() = cv[j].back() = ze;
                for (size_t y = 1; y < ny - 1; ++y)
                    cv[j][y] *= a0[j];

                for (size_t i = 0; i < nc; ++i) {
                    solver(ac[j][i], bc[j][i], cc[j][i], nv[i]);
                    std::transform(cv[j].begin(), cv[j].end(), nv[i].begin(), cv[j].begin(),
                                   [&c = aa[j][i]](const auto &a, const auto &b) { return a + c * b; });
                }

                _complete_mode_solution(j, x, nw, k0[j], ov, ph[j], cv[j]);
            };

            _compute(solve_func, nm, ny, nc, num_workers, buff_size);
        }

        Arg hx() const {
            return _hx;
        }

        Arg hy() const {
            return (_y1 - _y0) / (_ny - 1);
        }

        Arg hz() const {
            return _nz > 1 ? (_z1 - _z0) / (_nz - 1) : 0;
        }
        

    private:

        static constexpr auto im = Val(0, 1);
        static constexpr auto ze = Val(0);
        static constexpr auto on = Arg(1);
        static constexpr auto tw = Arg(2);

        BC _boundary_conditions;
        const size_t _nx{}, _ny{}, _nz{};
        coefficients<Val> _coefficients;
        const Arg _hx, _x0, _y0, _y1, _z0, _z1;

        static auto _get_thomas_solver(const size_t ny) {
            return [e=types::vector1d_t<Val>(ny), ny](
                    const types::vector1d_t<Val>& a,
                    const types::vector1d_t<Val>& b,
                    const types::vector1d_t<Val>& c,
                          types::vector1d_t<Val>& d) mutable {
                e[0] = c[0] / b[0];
                d[0] /= b[0];
                for (size_t i = 1; i < ny - 1; ++i) {
                    const auto w = b[i] - a[i] * e[i - 1];
                    e[i] = c[i] / w;
                    d[i] = (d[i] - a[i] * d[i - 1]) / w;
                }
                d.back() = (d.back() - a.back() * d[ny - 2]) / (b.back() - a.back() * e[ny - 2]);
                for (size_t i = ny - 1; i > 0; --i)
                    d[i - 1] -= e[i - 1] * d[i];
            };
        }

        static bool _all(const types::vector1d_t<bool>& values) {
            return std::all_of(values.begin(), values.end(), [](const auto& v) { return v; });
        }

        template<typename OV, typename PH, typename CV, typename K0>
        void _complete_mode_solution(const size_t& j, const Arg& x, const size_t& nw, const K0& k0, OV& ov, const PH& ph, const CV& cv) const {
            for (size_t i = 0, y = nw; i < _ny; ++i, ++y)
                for (size_t z = 0; z < _nz; ++z)
                    ov[i][z] = ph[i][z] * cv[y];

            on_mode_solution(j, x, k0, ov);

            for (size_t i = 0, y = nw; i < _ny; ++i, ++y) {
                const auto exp = std::exp(im * k0 * x);
                for (size_t z = 0; z < _nz; ++z)
                    ov[i][z] *= exp;
            }
        }

        template<typename SF>
        void _compute(const SF& solve_func, const size_t& nm, const size_t& ny, const size_t& nc, size_t num_workers, const size_t buff_size) const {
            threads::pool<size_t, Arg> pool(num_workers);

            types::vector3d_t<Val> ov(nm, types::vector2d_t<Val>(_ny, types::vector1d_t<Val>(_nz))),
                                   nv(nm, types::vector2d_t<Val>( nc, types::vector1d_t<Val>( ny)));

            std::vector solvers(nm, _get_thomas_solver(ny));
            threads::buffer_manager<types::vector2d_t<Val>> buffer(nm, buff_size, types::vector2d_t<Val>(_ny, types::vector1d_t<Val>(_nz, ze)));

            const auto& call = pool.add(
                [&, buffer_index=size_t(0)](const size_t& j, const size_t& i) {
                    while (true) {
                        auto data = buffer.current();

                        if (!data.ready())
                            break;

                        on_solution(_x0 + _hx * i, data.data());

                        for (size_t m = 0; m < _ny; ++m)
                            data.data()[m].assign(_nz, ze);

                        data.complete();
                        data.unlock();
                    }
                }
            );

            const auto& solution = pool.add(
                [&](const size_t& j, const size_t& i, const threads::task<size_t, Arg>& task) {
                    auto data = buffer.get(j);

                    if (!data.available()) {
                        task.push(j, i);
                        return;
                    }

                    data.unlock();

                    const auto x = _x0 + _hx * i;
                    solve_func(j, x, ov[j], nv[j], solvers[j]);

                    data.lock();
                    auto& ov_buff = data.data();

                    for (size_t m = 0; m < _ny; ++m)
                        for (size_t z = 0; z < _nz; ++z)
                            ov_buff[m][z] += ov[j][m][z];

                    if (data.complete())
                        call.push_single(j, i);

                    data.unlock();

                    if (i < _nx - 1)
                        task.push(j, i + 1);
                }
            );

            for (size_t j = 0; j < nm; ++j)
                solution.push(j, 1);

            pool.join();
        }

    };

}// namespace ample
