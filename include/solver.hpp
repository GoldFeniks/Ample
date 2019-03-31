#pragma once
#include <cmath>
#include <array>
#include <tuple>
#include <cstddef>
#include <complex>
#include <utility>
#include <cstring>
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/writer.hpp"
#include "utils/callback.hpp"

namespace acstc {

    using namespace std::complex_literals;

    template<typename Arg = types::real_t, typename Val = types::complex_t>
    class solver {

    public:

        solver(const Val& k0, const Arg& a,  const Arg& b,  const Arg& c,
               const Arg& x0, const Arg& x1, const Arg& nx,
               const Arg& y0, const Arg& y1, const Arg& ny) :
               _k0(k0), _a(a), _b(b), _c(c), _hx((x1 - x0) / (nx - 1)), _hy((y1 - y0) / (ny - 1)),
               _sq_k0(std::pow(k0, 2)), _sq_hy(std::pow(_hy, 2)), _nx(nx), _ny(ny) {
            const auto den = _sq_k0 * _sq_hy;
            const auto a0 = (_a - on) * im * _k0 * _hx;
            const auto a1 = (_b - _c) * im * _k0 * _hx;
            _g0 = tw + a0;
            _b0 = tw - a0;
            _g1 = tw * _c + a1;
            _b1 = tw * _c - a1;
            _ca = _g1 / den;
            _na = _b1 / den;
        }

        template<typename IV, typename KV, typename PV, typename FA>
        void operator()(const IV& init, const KV& ks, const PV& phi, FA&& callback_factory, const size_t k = 10, const size_t past_n = 0) const {
            types::vector1d_t<Val> nv(_ny), cv(_ny), ov(_ny);
            Val s0;
            types::vector1d_t<Val> va(_ny, _na), vb(_ny), cb(_ny), fv(_nx), lv(_nx), ss(_nx);

            for (const auto& pj : phi) {
                auto callback = callback_factory();
                for (const auto& [vv, kj] : utils::zip(init, ks)) {
                    const auto nb = std::pow(kj[0] / _k0, 2);

                    for (size_t m = 0; m < _ny; ++m) {
                        const auto dd = (std::pow(kj[m], 2) - _sq_k0 - tw / _sq_hy) / _sq_k0;
                        cb[m] = _g0 + _g1 * dd;
                        vb[m] = _b0 + _b1 * dd;
                    }

                    const auto de = on - _c * (on - nb);
                    const auto rr = tw * _k0 * _sq_hy / (_b - _c) / _hx;
                    const auto qq = tw * _c / _k0 / (_b - _c) / _hx;
                    const auto ka = _hx * _k0 * (_a - on - (_b - _c) * (on - nb)) / tw;
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
                    va[0] = va.back() = pq;
                    vb[0] = vb.back() = -ss[0];
                    s0 = -nq;

                    cv = vv;
                    for (size_t i = 0; i < _ny; ++i) {
                        ov[i] = cv[i] * pj[i];
                    }
                    callback(ov);

                    for (size_t n = 1; n < _nx; ++n) {
                        nv[0] = s0 * cv[1];
                        nv.back() = s0 * cv[_ny - 2];
                        for (size_t m = _start_index(n, past_n); m < n; ++m) {
                            nv[0] += fv[m] * ss[n - m];
                            nv.back() += lv[m] * ss[n - m];
                        }

                        for (size_t m = 1; m < _ny - 1; ++m)
                            nv[m] = (cv[m - 1] + cv[m + 1]) * _ca + cv[m] * cb[m];

                        _thomas_solver(va, vb, nv);
                        std::swap(cv, nv);

                        fv[n] = cv[0];
                        lv[n] = cv.back();

                        for (size_t i = 0; i < _ny; ++i)
                            ov[i] = cv[i] * pj[i];
                        callback(ov);
                    }
                }
            }
        }

    private:

        static constexpr auto im = Val(0, 1);
        static constexpr auto on = Arg(1);
        static constexpr auto tw = Arg(2);

        const Val _k0, _sq_k0;
        const Arg _a, _b, _c, _hx, _hy, _sq_hy;
        const size_t _nx, _ny;
        Val _g0, _g1, _b0, _b1, _ca, _na;

        static auto _start_index(const size_t n, const size_t m) {
            if (n < m)
                return size_t(1);
            return m ? n - m : size_t(1);
        }

        static auto _psqrt(const Val& value) {
            const auto res = std::sqrt(value);
            return res.real() < 0 ? -res : res;
        }

        auto _thomas_solver(const types::vector1d_t<Val>& a,
                const types::vector1d_t<Val>& b, types::vector1d_t<Val>& d) const {
            static types::vector1d_t<Val> c(_ny);
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
        }

    };

}// namespace acstc
