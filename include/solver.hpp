#pragma once
#include <cmath>
#include <array>
#include <tuple>
#include <mutex>
#include <thread>
#include <cstddef>
#include <complex>
#include <utility>
#include <cstring>
#include <condition_variable>
#include "config.hpp"
#include "io/writer.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/callback.hpp"
#include "utils/interpolation.hpp"

namespace acstc {
    
    using namespace std::complex_literals;

    template<typename Arg = types::real_t, typename Val = types::complex_t>
    class solver {

    public:

        solver(const Arg& a,  const Arg& b,  const Arg& c,
               const Arg& x0, const Arg& x1, const Arg& nx,
               const Arg& y0, const Arg& y1, const Arg& ny)
               : _a(a), _b(b), _c(c), _hx((x1 - x0) / (nx - 1)),
                 _hy((y1 - y0) / (ny - 1)), _sq_hy(std::pow(_hy, 2)), _x0(x0), _x1(x1), _y0(y0), _y1(y1), _nx(nx), _ny(ny) {}

        explicit solver(const config<Arg>& config) :
            solver(config.a(), config.b(), config.c(),
                   config.x0(), config.x1(), config.nx(),
                   config.y0(), config.y1(), config.ny()) {}

        template<typename IN, typename K0, typename CL>
        void solve(const IN& init,
                   const K0& k0,
                   const utils::linear_interpolated_data_2d<Arg, Val>& k_int,
                   const utils::linear_interpolated_data_2d<Arg, Arg>& phi_int,
                   CL&& callback,
                   const size_t past_n = 0,
                   const size_t num_workers = 1,
                   const size_t buff_size = 100) const {
            const auto mc = k0.size();
            if (init.size() != mc || k_int.size() != mc || phi_int.size() != mc)
                throw std::logic_error("Arguments init, k0, k, phi must be the same size");


            types::vector1d_t<Val> ca(mc), s0(mc), bv(_ny, Val(0)), g0(mc), g1(mc), b0(mc), b1(mc), sq_k0(mc);
            types::vector2d_t<Val> cv(mc, types::vector1d_t<Val>(_ny)),
                                   cb(mc, types::vector1d_t<Val>(_ny)),
                                   va(mc, types::vector1d_t<Val>(_ny)),
                                   vb(mc, types::vector1d_t<Val>(_ny)),
                                   fs(mc, types::vector1d_t<Val>(_nx)),
                                   ls(mc, types::vector1d_t<Val>(_nx)),
                                   fv(mc, types::vector1d_t<Val>(_nx)),
                                   lv(mc, types::vector1d_t<Val>(_nx));

            for (size_t j = 0; j < mc; ++j) {
                sq_k0[j] = std::pow(k0[j], 2);
                const auto den = sq_k0[j] * _sq_hy;
                const auto a0 = (_a - on) * im * k0[j] * _hx;
                const auto a1 = (_b - _c) * im * k0[j] * _hx;
                g0[j] = tw + a0;
                b0[j] = tw - a0;
                g1[j] = tw * _c + a1;
                b1[j] = tw * _c - a1;
                ca[j] = g1[j] / den;
                va[j].assign(_ny, b1[j] / den);
                s0[j] = im * tw * _c / k0[j] / (_b - _c) / _hx - on;
                va[j][0] = va[j].back() = s0[j] + tw;
                const auto k = k_int[j].template line(_x0, _y0, _y1, _ny);
                const auto phi = phi_int[j].template line(_x0, _y0, _y1, _ny);
                for (size_t i = 0; i < _ny; ++i) {
                    cv[j][i] = init[j][i];
                    bv[i] += phi[i] * init[j][i];
                }
                _fill_parameters(k0[j], k[0], fs[j], vb[j][0]);
                _fill_parameters(k0[j], k.back(), ls[j], vb[j].back());
            }

            callback(bv);

            auto solve_func = [&](const size_t j0, const size_t j1, auto&& call) {
                types::vector1d_t<Val> nv(_ny), ov(_ny, Val(0)), k(_ny);
                types::vector1d_t<Arg> phi(_ny);

                auto x = _hx;
                auto solver = _get_thomas_solver();

                for (size_t i = 1; i < _nx; ++i) {
                    ov.assign(_ny, Val(0));
                    for (size_t j = j0; j < j1; ++j) {
                        k_int[j].template line(x, _y0, _y1, k);
                        phi_int[j].template line(x, _y0, _y1, phi);

                        for (size_t m = 1; m < _ny - 1; ++m) {
                            const auto dd = (std::pow(k[m], 2) - sq_k0[j] - tw / _sq_hy) / sq_k0[j];
                            cb[j][m] = g0[j] + g1[j] * dd;
                            vb[j][m] = b0[j] + b1[j] * dd;
                        }

                        nv[0] = s0[j] * cv[j][1];
                        nv.back() = s0[j] * cv[j][_ny - 2];
                        for (size_t m = _start_index(i, past_n); m < i; ++m) {
                            nv[0] += fv[j][m] * fs[j][i - m];
                            nv.back() += lv[j][m] * ls[j][i - m];
                        }

                        for (size_t m = 1; m < _ny - 1; ++m)
                            nv[m] = (cv[j][m - 1] + cv[j][m + 1]) * ca[j] + cv[j][m] * cb[j][m];

                        solver(va[j], vb[j], nv);
                        std::swap(cv[j], nv);

                        fv[j][i] = cv[j][0];
                        lv[j][i] = cv[j].back();

                        for (size_t m = 0; m < _ny; ++m)
                            ov[m] += phi[m] * cv[j][m] * std::exp(im * k0[j] * x);
                    }

                    call(ov);
                    x += _hx;
                }
            };

            _compute(solve_func, callback, mc, num_workers, buff_size);

        }

        template<typename IN, typename K0, typename CL>
        void solve(const IN& init,
                         const K0& k0,
                         const utils::linear_interpolated_data_1d<Arg, Val>& k_int,
                         const utils::linear_interpolated_data_1d<Arg, Arg>& phi_int,
                         CL&& callback,
                         const size_t past_n = 0,
                         const size_t num_workers = 1,
                         const size_t buff_size = 100) const {
            const auto mc = k0.size();
            if (init.size() != mc || k_int.size() != mc || phi_int.size() != mc)
                throw std::logic_error("Vectors init, k0, k, phi must be the same size");

            types::vector2d_t<Val> k(mc);
            types::vector2d_t<Arg> phi(mc);
            for (size_t j = 0; j < mc; ++j) {
                k[j] = k_int[j].line(_y0, _y1, _ny);
                phi[j] = phi_int[j].line(_y0, _y1, _ny);
            }

            types::vector1d_t<Val> ca(mc), s0(mc), bv(_ny, Val(0));
            types::vector2d_t<Val> cv(mc, types::vector1d_t<Val>(_ny)),
                                   cb(mc, types::vector1d_t<Val>(_ny)),
                                   va(mc, types::vector1d_t<Val>(_ny)),
                                   vb(mc, types::vector1d_t<Val>(_ny)),
                                   fs(mc, types::vector1d_t<Val>(_nx)),
                                   ls(mc, types::vector1d_t<Val>(_nx)),
                                   fv(mc, types::vector1d_t<Val>(_nx)),
                                   lv(mc, types::vector1d_t<Val>(_nx));

            for (size_t j = 0; j < mc; ++j) {
                const auto sq_k0 = std::pow(k0[j], 2);
                const auto den = sq_k0 * _sq_hy;
                const auto a0 = (_a - on) * im * k0[j] * _hx;
                const auto a1 = (_b - _c) * im * k0[j] * _hx;
                const auto g0 = tw + a0;
                const auto b0 = tw - a0;
                const auto g1 = tw * _c + a1;
                const auto b1 = tw * _c - a1;
                ca[j] = g1 / den;
                va[j].assign(_ny, b1 / den);
                s0[j] = im * tw * _c / k0[j] / (_b - _c) / _hx - on;
                va[j][0] = va[j].back() = s0[j] + tw;
                for (size_t i = 0; i < _ny; ++i) {
                    const auto dd = (std::pow(k[j][i], 2) - sq_k0 - tw / _sq_hy) / sq_k0;
                    cb[j][i] = g0 + g1 * dd;
                    vb[j][i] = b0 + b1 * dd;
                    cv[j][i] = init[j][i];
                    bv[i] += phi[j][i] * init[j][i];
                }
                _fill_parameters(k0[j], k[j][0], fs[j], vb[j][0]);
                _fill_parameters(k0[j], k[j].back(), ls[j], vb[j].back());
            }

            callback(bv);

            auto solve_func = [&](const size_t j0, const size_t j1, auto&& call) {
                types::vector1d_t<Val> nv(_ny), ov(_ny, Val(0));

                auto x = _hx;
                auto solver = _get_thomas_solver();

                for (size_t i = 1; i < _nx; ++i) {
                    ov.assign(_ny, Val(0));
                    for (size_t j = j0; j < j1; ++j) {
                        nv[0] = s0[j] * cv[j][1];
                        nv.back() = s0[j] * cv[j][_ny - 2];
                        for (size_t m = _start_index(i, past_n); m < i; ++m) {
                            nv[0] += fv[j][m] * fs[j][i - m];
                            nv.back() += lv[j][m] * ls[j][i - m];
                        }

                        for (size_t m = 1; m < _ny - 1; ++m)
                            nv[m] = (cv[j][m - 1] + cv[j][m + 1]) * ca[j] + cv[j][m] * cb[j][m];

                        solver(va[j], vb[j], nv);
                        std::swap(cv[j], nv);

                        fv[j][i] = cv[j][0];
                        lv[j][i] = cv[j].back();

                        for (size_t m = 0; m < _ny; ++m)
                            ov[m] += phi[j][m] * cv[j][m] * std::exp(im * k0[j] * x);
                    }
                    call(ov);
                    x += _hx;
                }
            };

            _compute(solve_func, callback, mc, num_workers, buff_size);
        }

    private:

        static constexpr auto im = Val(0, 1);
        static constexpr auto on = Arg(1);
        static constexpr auto tw = Arg(2);

        const Arg _a, _b, _c, _hx, _hy, _sq_hy, _x0, _x1, _y0, _y1;
        const size_t _nx, _ny;

        static auto _start_index(const size_t n, const size_t m) {
            if (n < m)
                return size_t(1);
            return m ? n - m : size_t(1);
        }

        static auto _psqrt(const Val& value) {
            const auto res = std::sqrt(value);
            return res.real() < 0 ? -res : res;
        }

        void _fill_parameters(const Val& k0, const Val& k, types::vector1d_t<Val>& ss, Val& vb) const {
            const auto nb = std::pow(k / k0, 2);
            const auto de = on - _c * (on - nb);
            const auto rr = tw * k0 * _sq_hy / (_b - _c) / _hx;
            const auto qq = tw * _c / k0 / (_b - _c) / _hx;
            const auto ka = _hx * k0 * (_a - on - (_b - _c) * (on - nb)) / tw;
            const auto ga = rr * de;
            const auto si = ka * -rr;
            const auto pg = ga + im * si;
            const auto ng = ga - im * si;
            const auto gq = ga - Arg(4) * qq;
            const auto s4 = si + Arg(4);
            const auto aa = pg * (gq + im * s4);
            const auto bb = ga * gq + si * s4;
            const auto cc = ng * (gq - im * s4);
            const auto sa = _psqrt(aa);
            const auto sc = _psqrt(cc);
            const auto la = sa / sc;
            const auto mu = bb / sa / sc;
            const auto pq = on + im * qq;
            const auto nq = on - im * qq;
            const auto xi = std::arg(nq / pq);
            const auto ph = std::arg(aa);
            const auto bh = pg + ng * std::exp(-im * xi);
            const auto ze = im * std::sqrt(std::abs(aa)) * std::exp(im * ph / tw) / tw;
            const auto be = nq + im * ng / tw + cc * (pq - im * pg / tw) / bb;

            ss[0] = pq - im * (pg + sa) / tw;
            ss[1] = nq + im * (ng + bb / sa) / tw;
            ss[2] = mu * (ss[1] + ss[0] / la / mu - be) / la / tw;
            for (size_t n = 2; n < _nx - 1; ++n)
                ss[n + 1] = (Arg(2 * n - 1) * mu * ss[n] / la -
                                Arg(n - 2) * ss[n - 1] / std::pow(la, 2)) / Arg(n + 1);
            vb = -ss[0];
        }

        auto _get_thomas_solver() const {
            return [this, c=types::vector1d_t<Val>(_ny)](const types::vector1d_t<Val>& a,
                    const types::vector1d_t<Val>& b, types::vector1d_t<Val>& d) mutable {
                c[0] = a[0] / b[0];
                d[0] /= b[0];
                for (size_t i = 1; i < _ny - 1; ++i) {
                    const auto w = b[i] - a[i] * c[i - 1];
                    c[i] = a[i] / w;
                    d[i] = (d[i] - a[i] * d[i - 1]) / w;
                }
                d.back() = (d.back() - a.back() * d[_ny - 2]) / (b.back() - a.back() * c[_ny - 2]);
                for (size_t i = _ny - 1; i > 0; --i)
                    d[i - 1] -= c[i - 1] * d[i];
            };
        }

        template<typename SF, typename CL>
        void _compute(const SF& solve_func, CL&& callback, const size_t mc, size_t num_workers, const size_t buff_size) const {
            num_workers = std::min(mc, num_workers);

            if (num_workers <= 1) {
                solve_func(0, mc, callback);
                return;
            }

            types::vector1d_t<std::thread> workers;
            workers.reserve(num_workers);
            types::vector2d_t<Val> ov_buff(buff_size, types::vector1d_t<Val>(_ny, Val(0)));
            types::vector2d_t<bool> done_buff(buff_size, types::vector1d_t<bool>(num_workers, false));
            types::vector1d_t<std::mutex> buff_mutex(buff_size);
            const auto mpw = mc / num_workers;

            for (size_t i = 0; i < num_workers; ++i)
                workers.emplace_back([&, i](){
                    solve_func(mpw * i, i == num_workers - 1 ? mc : mpw * (i + 1), [&, i, in=size_t(0)](const auto& data) mutable {
                        while (true) {
                            buff_mutex[in].lock();
                            if (!done_buff[in][i])
                                break;
                            std::cout << "Buffer size exceeded" << std::endl;
                            buff_mutex[in].unlock();
                        }
                        for (size_t m = 0; m < _ny; ++m)
                            ov_buff[in][m] += data[m];
                        done_buff[in][i] = true;
                        buff_mutex[in].unlock();
                        in = (in + 1) % buff_size;
                    });
                });

            size_t in = 1, bi = 0;
            while (in < _nx) {
                buff_mutex[bi].lock();
                if (std::all_of(done_buff[bi].begin(), done_buff[bi].end(), [](const auto& val) { return val; })) {
                    callback(ov_buff[bi]);
                    done_buff[bi].assign(num_workers, false);
                    ov_buff[bi].assign(_ny, Val(0));
                    bi = (bi + 1) % buff_size;
                    ++in;
                }
                buff_mutex[bi].unlock();
            }

            for (auto& it : workers)
                it.join();
        }

    };

}// namespace acstc
