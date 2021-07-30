#pragma once
#include <cmath>
#include <mutex>
#include <thread>
#include <cstddef>
#include <complex>
#include <algorithm>
#include <condition_variable>
#include "config.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/assert.hpp"
#include "coefficients.hpp"
#include "boundary_conditions.hpp"
#include "utils/interpolation.hpp"

//WINDOWS DEFINES VERY NICE MACRO FOR NO REASON
#ifdef IN
#undef IN
#endif

namespace ample {
    
    using namespace std::complex_literals;

    template<typename BC, typename Arg = typename BC::arg_t, typename Val = typename BC::val_t>
    class solver {

    public:

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

        template<typename IN, typename CL, typename VL>
        void solve(const IN& init,
                   const types::vector1d_t<VL>& k0,
                   const utils::linear_interpolated_data_1d<Arg, VL>& k_int,
                   const utils::linear_interpolated_data_2d<Arg, Arg>& phi_int,
                   CL&& callback,
                   const size_t num_workers = 1,
                   const size_t buff_size = 100) const {
            const auto nm = k0.size();
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

            callback(_x0, bv);

            auto solve_func = [&, &ac=ac, &bc=bc, &cc=cc](const size_t j0, const size_t j1, auto&& call) {
                types::vector2d_t<Val> ov(_ny, types::vector1d_t<Val>(_nz)),
                                       nv( nc, types::vector1d_t<Val>( ny));

                auto x = _x0 + _hx;
                auto solver = _get_thomas_solver(ny);

                for (size_t _ = 1; _ < _nx; ++_) {
                    for (size_t y = 0; y < _ny; ++y)
                        ov[y].assign(_nz, ze);

                    for (size_t j = j0; j < j1; ++j) {
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

                        for (size_t i = 0, y = nw; i < _ny; ++i, ++y) {
                            const auto exp = cv[j][y] * std::exp(im * k0[j] * x);
                            for (size_t z = 0; z < _nz; ++z)
                                ov[i][z] += ph[j][i][z] * exp;
                        }
                    }

                    call(x, ov);
                    x += _hx;
                }
            };

            _compute(solve_func, callback, nm, num_workers, buff_size);
        }

        template<typename IN, typename CL, typename VL>
        void solve(const IN& init,
                   const types::vector1d_t<VL>& k0,
                   const utils::linear_interpolated_data_2d<Arg, VL>& k_int,
                   const utils::linear_interpolated_data_3d<Arg, Arg>& phi_int,
                   CL&& callback,
                   const size_t num_workers = 1,
                   const size_t buff_size = 100) const {
            const auto nm = k0.size();
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

            callback(_x0, bv);

            types::vector2d_t<VL>  kk(nm, types::vector1d_t<VL> (_ny));
            types::vector3d_t<Arg> ph(nm, types::vector2d_t<Arg>(_ny, types::vector1d_t<Arg>(_nz)));

            auto solve_func = [&, &ac=ac, &bc=bc, &cc=cc](const size_t j0, const size_t j1, auto&& call) {
                types::vector2d_t<Val> ov(_ny, types::vector1d_t<Val>(_nz)),
                                       nv(nc, types::vector1d_t<Val>( ny));

                auto x = _x0 + _hx;
                auto solver = _get_thomas_solver(ny);

                for (size_t _ = 1; _ < _nx; ++_) {
                    for (size_t y = 0; y < _ny; ++y)
                        ov[y].assign(_nz, ze);

                    for (size_t j = j0; j < j1; ++j) {
                        k_int[j].line(x, _y0, _y1, kk[j]);
                        phi_int[j].field(x, _y0, _y1, _z0, _z1, ph[j]);
                    }

                    band_builder.update(kk, j0, j1);

                    for (size_t j = j0; j < j1; ++j) {
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

                        for (size_t i = 0, y = nw; i < _ny; ++i, ++y) {
                            const auto exp = cv[j][y] * std::exp(im * k0[j] * x);
                            for (size_t z = 0; z < _nz; ++z)
                                ov[i][z] += ph[j][i][z] * exp;
                        }
                    }

                    call(x, ov);
                    x += _hx;
                }
            };

            _compute(solve_func, callback, nm, num_workers, buff_size);
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

        template<typename SF, typename CL>
        void _compute(const SF& solve_func, CL&& callback, const size_t mc, size_t num_workers, const size_t buff_size) const {
            num_workers = std::min(mc + 1, num_workers);

            if (num_workers <= 1) {
                solve_func(0, mc, callback);
                return;
            }

            --num_workers; // one worker is used to output data

            types::vector1d_t<std::mutex> buff_mutex(buff_size);
            types::vector1d_t<std::condition_variable> buff_cv(buff_size);
            types::vector2d_t<bool> done_buff(buff_size, types::vector1d_t<bool>(num_workers, false));
            types::vector3d_t<Val> ov_buff(buff_size, types::vector2d_t<Val>(_ny, types::vector1d_t<Val>(_nz, ze)));

            const auto mpw = mc / num_workers;

            types::vector1d_t<std::thread> workers;
            workers.reserve(num_workers);

            for (size_t i = 0; i < num_workers; ++i)
                workers.emplace_back([&, i](){
                    solve_func(mpw * i, i == num_workers - 1 ? mc : mpw * (i + 1), [&, i, in=size_t(0)](const auto& x, const auto& data) mutable {
                        std::unique_lock<std::mutex> lk(buff_mutex[in]);
                        if (done_buff[in][i])
                            buff_cv[in].wait(lk, [&]{ return !done_buff[in][i]; });

                        for (size_t m = 0; m < _ny; ++m)
                            for (size_t z = 0; z < _nz; ++z)
                                ov_buff[in][m][z] += data[m][z];

                        done_buff[in][i] = true;
                        if (_all(done_buff[in])) {
                            lk.unlock();
                            buff_cv[in].notify_all();
                        } else
                            lk.unlock();

                        in = (in + 1) % buff_size;
                    });
                });

            for (size_t in = 1, bi = 0; in < _nx; ++in, bi = (bi + 1) % buff_size) {
                std::unique_lock<std::mutex> lk(buff_mutex[bi]);
                if (!_all(done_buff[bi]))
                    buff_cv[bi].wait(lk, [&]{ return _all(done_buff[bi]); });

                callback(_x0 + in * _hx, ov_buff[bi]);

                for (size_t y = 0; y < _ny; ++y)
                    ov_buff[bi][y].assign(_nz, ze);

                done_buff[bi].assign(num_workers, false);
                lk.unlock();
                buff_cv[bi].notify_all();
            }

            for (auto& it : workers)
                it.join();
        }

    };

}// namespace ample
