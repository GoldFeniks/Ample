#pragma once
#include <cmath>
#include <array>
#include <tuple>
#include <cstddef>
#include <complex>
#include <utility>
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
            const auto [b0, b1, g0, g1, ca, na] = _calc_coefficients(coefficients, k0, sq_k0, hx, hy, sq_hy);
            return [sq_k0, sq_hy, b0, b1, g0, g1, ca, na, nx, ny]
                (const types::vector2d_t<Val>& init, const types::vector1d_t<Arg>& k, utils::basic_writer<Val>& writer) {
                types::vector1d_t<Val> nv(ny);
                auto cb = ca;
                auto nb = na;

                for (size_t j = 0; j < init.size(); ++j) {
                    const auto d = (Arg(2) / sq_hy + sq_k0 - std::pow(k[j], 2)) / sq_k0;
                    for (size_t i = 0; i < b0.size(); ++i) {
                        cb[i] = g0[i] + g1[i] * d;
                        nb[i] = b0[i] + g0[i] * d;
                    }

                    auto cv = init[j];
                    writer.write(cv);
                    for (size_t n = 1; n < nx; ++n) {
                        for (size_t i = 0; i < b0.size(); ++i) {
                            const auto a = ca[i];
                            const auto c = cb[i];
                            nv[0] = cv[0] * c + cv[1] * a;
                            nv.back() = cv[ny - 2] * a + cv.back() * c;
                            for (size_t m = 1; m < ny - 1; ++m)
                                nv[m] = a * (cv[m - 1] + cv[m + 1]) + cv[m] * c;
                            _thomas_solve(na[i], nb[i], nv);
                            std::swap(cv, nv);
                        }
                        writer.write(cv);
                    }
                }
            };
        }

    private:

        template<typename A>
        static auto _calc_coefficients(const A& coefficients, const Arg& k0, const Arg& sq_k0,
                                       const Arg& hx, const Arg& hy, const Arg& sq_hy) {
            std::array<Val, A::size> g0, g1, b0, b1, ca, na;
            const auto den = sq_k0 * sq_hy;
            const auto im = Val(1i);
            for (size_t i = 0; i < A::size; ++i) {
                const auto& [a, b, c] = coefficients[i];
                const auto a0 = a - Arg(1);
                const auto a1 = c - b;
                g0[i] = Arg(2) + im * k0 * hx * a0;
                b0[i] = Arg(2) - im * k0 * hx * a0;
                g1[i] = Arg(2) * c + im * k0 * hx * a1;
                b1[i] = Arg(2) * c - im * k0 * hx * a1;
                ca[i] = -g1[i] / den;
                na[i] = -b1[i] / den;
            }
            return std::tuple(g0, g1, b0, b1, ca, na);
        }

        static void _thomas_solve(const Val& a, const Val& b, types::vector1d_t<Val>& f) {
            types::vector1d_t<Val> buff(f.size());
            buff[0] = b;
            for (size_t i = 1; i < f.size(); ++i) {
                const auto w = a / buff[i - 1];
                buff[i] -= w * a;
                f[i] -= w * f[i - 1];
            }
            f.back() /= buff.back();
            auto i = f.size() - 1;
            do {
                --i;
                f[i] = (f[i] - a * f[i + 1]) / buff[i];
            } while (i);
        }
    };

}// namespace acstc
