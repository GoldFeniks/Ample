#pragma once
#include <cmath>
#include <array>
#include <tuple>
#include <cstddef>
#include <complex>
#include <utility>
#include <cstring>
#include "utils/types.hpp"
#include "utils/writer.hpp"

namespace acstc {

    using namespace std::complex_literals;

    template<typename Arg = types::real_t, typename Val = types::complex_t>
    class solver {

    public:

        solver() = delete;

        template<typename A>
        static auto create(const Arg& k0, const A& coefficients,
                           const Arg& x0, const Arg& x1, const Arg& nx,
                           const Arg& y0, const Arg& y1, const Arg& ny) {
            const auto hx = (x1 - x0) / (nx - 1);
            const auto hy = (y1 - y0) / (ny - 1);
            const auto sq_k0 = std::pow(k0, 2);
            const auto sq_hy = std::pow(hy, 2);
            const auto [g0, g1, b0, b1, ca, na] = _calc_coefficients(coefficients, k0, sq_k0, hx, hy, sq_hy);
            return [nx, ny, k0, sq_k0, hx, hy, sq_hy, b0, b1, g0, g1, ca, na, coefficients]
                (const auto& init, const auto& k, auto&& callback, const size_t past_n = 0) {
                types::vector1d_t<Val> nv(ny);
                std::array<Val, A::size> cb, s0;
                std::array<types::vector1d_t<Val>, A::size> ss, va, vb, fv, lv;

                for (size_t i = 0; i < va.size(); ++i) {
                    va[i].resize(ny, na[i]);
                    vb[i].resize(ny);
                    ss[i].resize(nx);
                    fv[i].resize(nx);
                    lv[i].resize(nx);
                }

                auto thomas_solver = _get_thomas_solver(ny);

                for (size_t j = 0; j < init.size(); ++j) {
                    const auto nb = k[j] / k0;
                    const auto dd = (std::pow(k[j], 2) - sq_k0 - tw / sq_hy) / sq_k0;
                    for (size_t i = 0; i < b0.size(); ++i) {
                        const auto& [a, b, c] = coefficients[i];

                        cb[i] = g0[i] + g1[i] * dd;
                        vb[i].assign(ny, b0[i] + b1[i] * dd);

                        const auto de = on - c * (on - nb);
                        const auto rr = tw * k0 * sq_hy / (b - c) / hx;
                        const auto qq = tw * c / k0 / (b - c) / hx;
                        const auto ka = hx * k0 * (a - on - (b - c) * (on - nb)) / tw;
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

                        ss[i][0] = pq - im * (pg + sa) / tw;
                        ss[i][1] = nq + im * (ng + bb / sa) / tw;
                        ss[i][2] = mu * (ss[i][1] + ss[i][0] / la / mu - be) / la / tw;
                        for (size_t n = 2; n < nx - 1; ++n)
                            ss[i][n + 1] = (Arg(2 * n - 1) * mu * ss[i][n] / la -
                                            Arg(n - 2) * ss[i][n - 1] / std::pow(la, 2)) / Arg(n + 1);
                        va[i][0] = va[i].back() = pq;
                        vb[i][0] = vb[i].back() = -ss[i][0];
                        s0[i] = -nq;
                    }

                    auto cv = init[j];
                    callback(cv);

                    for (size_t n = 1; n < nx; ++n) {
                        for (size_t i = 0; i < b0.size(); ++i) {
                            nv[0] = s0[i] * cv[1];
                            nv.back() = s0[i] * cv[ny - 2];
                            for (size_t m = _start_index(n, past_n); m < n; ++m) {
                                nv[0] += fv[i][m] * ss[i][n - m];
                                nv.back() += lv[i][m] * ss[i][n - m];
                            }

                            const auto a = ca[i];
                            const auto b = cb[i];
                            for (size_t m = 1; m < ny - 1; ++m)
                                nv[m] = (cv[m - 1] + cv[m + 1]) * a + cv[m] * b;

                            thomas_solver(va[i], vb[i], nv);
                            std::swap(cv, nv);

                            fv[i][n] = cv[0];
                            lv[i][n] = cv.back();
                        }
                        callback(cv);
                    }
                }
            };
        }

    private:

        static constexpr auto im = Val(0, 1);
        static constexpr auto on = Arg(1);
        static constexpr auto tw = Arg(2);

        static auto _start_index(const size_t n, const size_t m) {
            if (n < m)
                return size_t(1);
            return m ? n - m : size_t(1);
        }

        static auto _psqrt(const Val& value) {
            const auto res = std::sqrt(value);
            return res.real() < 0 ? -res : res;
        }

        template<typename A>
        static auto _calc_coefficients(const A& coefficients, const Arg& k0, const Arg& sq_k0,
                                       const Arg& hx, const Arg& hy, const Arg& sq_hy) {
            std::array<Val, A::size> g0, g1, b0, b1, ca, na;
            const auto den = sq_k0 * sq_hy;
            for (size_t i = 0; i < A::size; ++i) {
                const auto& [a, b, c] = coefficients[i];
                const auto a0 = (a - on) * im * k0 * hx;
                const auto a1 = (b - c) * im * k0 * hx;
                g0[i] = tw + a0;
                b0[i] = tw - a0;
                g1[i] = tw * c + a1;
                b1[i] = tw * c - a1;
                ca[i] = g1[i] / den;
                na[i] = b1[i] / den;

            }
            return std::make_tuple(g0, g1, b0, b1, ca, na);
        }

        static auto _get_thomas_solver(const size_t n) {
            return [n, c=types::vector1d_t<Val>(n)]
                (const types::vector1d_t<Val>& a, const types::vector1d_t<Val>& b, types::vector1d_t<Val>& d) mutable {
                    c[0] = a[0] / b[0];
                    d[0] /= b[0];
                    for (size_t i = 1; i < n - 1; ++i) {
                        const auto w = b[i] - a[i] * c[i - 1];
                        c[i] = a[i] / w;
                        d[i] = (d[i] - a[i] * d[i - 1]) / w;
                    }
                    d.back() = (d.back() - a.back() * d[n - 2]) / (b.back() - a.back() * c[n - 2]);
                    for (size_t i = n - 1; i > 0; --i)
                        d[i - 1] -= c[i - 1] * d[i];
            };
        }

    };

}// namespace acstc
