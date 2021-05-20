#pragma once

#include <cmath>
#include <tuple>
#include <complex>
#include <cstddef>
#include <algorithm>
#include <functional>
#include <type_traits>
#include "feniks/zip.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/assert.hpp"
#include "utils/interpolation.hpp"

namespace ample {

    namespace _impl {

        template<typename T, typename V>
        auto find_less(const T& value, const V& values, size_t i) {
            while (i > 1 && value < values[i - 1])
                --i;
            return i - 1;
        }

        template<typename T, typename V>
        auto find_greater(const T& value, const V& values, size_t i) {
            const auto m = values.size() - 2;
            while (i < m && value > values[i + 1])
                ++i;
            return i;
        }

        template<typename T, typename V>
        inline auto diff(const size_t& i, const T& h, const V& data) {
            if (i == 0)
                return (-T(3) * data[0] + T(4) * data[1] - data[2]) / (T(2) * h);
            if (i == data.size() - 1)
                return (data[i - 2] - T(4) * data[i - 1] + T(3) * data[i]) / (T(2) * h);
            return (data[i + 1] - data[i - 1]) / (T(2) * h);
        }

        template<typename T, typename C, typename V>
        void _taper_interval(const size_t& l, const size_t& r, const T& rc, const T& d, const C& coords, V& values) {
            for (size_t i = l + 1; i < r; ++i) {
                const auto v = std::pow(rc - coords[i], 2) / d;
                values[i] *= std::exp(v * v * (T(1) / (v - T(1)) - T(1)));
            }
        }

        template<typename C, typename V>
        void taper_interval(const size_t& il, const size_t& jl, const size_t& ir, const size_t& jr, const C& coords, V& values) {
            const auto lc = coords[jl];
            const auto rc = coords[jr];

            _taper_interval(il, jl, lc, std::pow(lc - coords[il], 2), coords, values);
            _taper_interval(jr, ir, rc, std::pow(rc - coords[ir], 2), coords, values);

            constexpr std::remove_reference_t<decltype(values[0])> zero(0);

            values[il] = values[ir] = zero;
        }

        template<typename C, typename V>
        void taper_interval_all(const size_t& il, const size_t& jl, const size_t& ir, const size_t& jr, const C& coords, V& values) {
            const auto lc = coords[jl];
            const auto rc = coords[jr];

            constexpr std::remove_reference_t<decltype(values[0][0])> zero(0);

            for (auto& it : values) {
                _taper_interval(il, jl, lc, std::pow(lc - coords[il], 2), coords, it);
                _taper_interval(jr, ir, rc, std::pow(rc - coords[ir], 2), coords, it);

                it[il] = it[ir] = zero;
            }
        }

        template<typename T, typename C, typename V>
        void _taper_tail(const size_t& l, const size_t& r, const T& rc, const T& d, const C& coords, V& values) {
            for (size_t i = l; i < r; ++i)
                values[i] *= std::exp(-std::pow(rc - coords[i], 2) / d);
        }

        template<typename C, typename V>
        void taper_tail(const size_t& il, const size_t& jl, const size_t& ir, const size_t& jr, const C& coords, V& values) {
            const auto lc = coords[il];
            const auto rc = coords[ir];

            _taper_tail(0,                 il, lc, std::pow(lc - coords[jl], 2), coords, values);
            _taper_tail(ir + 1, coords.size(), rc, std::pow(rc - coords[jr], 2), coords, values);
        }

        template<typename C, typename V>
        void taper_tail_all(const size_t& il, const size_t& jl, const size_t& ir, const size_t& jr, const C& coords, V& values) {
            const auto lc = coords[il];
            const auto rc = coords[ir];

            const auto dl = std::pow(lc - coords[jl], 2);
            const auto dr = std::pow(rc - coords[jr], 2);

            for (auto& it : values) {
                _taper_tail(0,                 il, lc, dl, coords, it);
                _taper_tail(ir + 1, coords.size(), rc, dr, coords, it);
            }
        }

    }// namespace _impl

    template<typename T>
    class percentage_tapering {

    public:

        explicit percentage_tapering(const T& p) : _pl(p), _pr(p) {}
        explicit percentage_tapering(const T& pl, const T& pr) : _pl(pl), _pr(pr) {}

        template<typename C, typename A, typename V>
        void interval(const size_t& il, const size_t& ir, const C& coords, const A&, V& values) const {
            const auto [jl, jr] = _get_boundaries(il, ir);

            _impl::taper_interval(il, jl, ir, jr, coords, values);
        }

        template<typename C, typename A, typename V>
        void interval_all(const size_t& il, const size_t& ir, const C& coords, const A&, V& values) const {
            const auto [jl, jr] = _get_boundaries(il, ir);

            _impl::taper_interval_all(il, jl, ir, jr, coords, values);
        }

        template<typename C, typename A, typename V>
        void tail(const size_t& il, const size_t& ir, const C& coords, const A&, V& values) const {
            const auto [jl, jr] = _get_boundaries(il, ir);

            _impl::taper_tail(il, jl, ir, jr, coords, values);
        }

        template<typename C, typename A, typename V>
        void tail_all(const size_t& il, const size_t& ir, const C& coords, const A&, V& values) const {
            const auto [jl, jr] = _get_boundaries(il, ir);

            _impl::taper_tail_all(il, jl, ir, jr, coords, values);
        }

    private:

        const T _pl, _pr;

        auto _get_boundaries(const size_t& il, const size_t& ir) const {
            const auto id = ir - il;
            return std::make_tuple(il + size_t(id * _pl), ir - size_t(id * _pr));
        }

    };

    template<typename T>
    class angled_tapering {

    public:

        explicit angled_tapering(const T& w) : _wl(w), _wr(w) {}
        angled_tapering(const T& wl, const T& wr) : _wl(wl), _wr(wr) {}

        template<typename C, typename A, typename V>
        void interval(const size_t& il, const size_t& ir, const C&, const A& angles, V& values) const {
            const auto [jl, jr] = _get_boundaries(il, ir, angles);

            _impl::taper_interval(il, jl, ir, jr, angles, values);
        }

        template<typename C, typename A, typename V>
        void interval_all(const size_t& il, const size_t& ir, const C&, const A& angles, V& values) const {
            const auto [jl, jr] = _get_boundaries(il, ir, angles);

            _impl::taper_interval_all(il, jl, ir, jr, angles, values);
        }

        template<typename C, typename A, typename V>
        void tail(const size_t& il, const size_t& ir, const C&, const A& angles, V& values) const {
            const auto [jl, jr] = _get_boundaries(il, ir, angles);

            _impl::taper_tail(il, jl, ir, jr, angles, values);
        }

        template<typename C, typename A, typename V>
        void tail_all(const size_t& il, const size_t& ir, const C&, const A& angles, V& values) const {
            const auto [jl, jr] = _get_boundaries(il, ir, angles);

            _impl::taper_tail_all(il, jl, ir, jr, angles, values);
        }

    private:

        const T _wl, _wr;

        template<typename A>
        auto _get_boundaries(const size_t& il, const size_t& ir, const A& angles) const {
            return std::make_tuple(
                std::distance(angles.begin(), std::upper_bound(angles.begin(), angles.end(), angles[il] + _wl)),
                std::distance(angles.begin(), std::lower_bound(angles.begin(), angles.end(), angles[ir] - _wr))
            );
        }

    };

    using namespace std::placeholders;

    template<typename Arg, typename Val = Arg>
    class initial_conditions {

    public:

        using func_t = std::function<types::vector2d_t<Val>(const Arg&, const Arg&, const size_t&)>;

        explicit initial_conditions(func_t function) : _function(std::move(function)) {}

        template<typename F, typename... Args>
        static initial_conditions create(F func, Args... args) {
            return initial_conditions(std::bind(&_create<F, Args...>, _1, _2, _3, std::move(func), std::move(args)...));
        }

        types::vector2d_t<Val> make(const Arg& a, const Arg& b, const size_t& n, const size_t& m) const {
            auto result = _function(a, b, n);
            utils::dynamic_assert(result.size() >= m,
                "Initial conditions function returned vector of insufficient size(", result.size(), "), expected at least ", m);
            for (size_t i = 0; i < m; ++i)
                utils::dynamic_assert(result[i].size() == n,
                    "Initial conditions function return vector of incorrect size(", result[i].size(), ") at level ", i, ", expected size ", n);

            return result;
        }

    private:

        func_t _function;

        template<typename F, typename... Args>
        static types::vector2d_t<Val> _create(const Arg& a, const Arg& b, const size_t& n, const F& func, const Args&... args) {
            const auto xs = utils::mesh_1d(a, b, n);

            const auto zip = feniks::zip(std::forward<const Args>(args)...);
            types::vector2d_t<Val> result(zip.size());
            for (size_t i = 0; i < result.size(); ++i) {
                result[i].reserve(xs.size());
                auto tuple = std::tuple_cat(std::make_tuple(xs[0]), zip[i]);
                for (size_t j = 0; j < xs.size(); ++j) {
                    std::get<0>(tuple) = xs[j];
                    result[i].emplace_back(std::apply(func, tuple));
                }
            }
            return result;
        }

    };

    template<typename Val, typename Arg, typename AV, typename WV>
    auto generalized_gaussian_source(const Arg& a, const Arg& b, const size_t& n, const Arg& x0,
            const Arg& t1, const Arg& t2, const AV& ca, const WV& cw) {
        return initial_conditions<Arg, Val>::create(
            [x0, t1, t2, im=Val(0, 1)](Arg x, const auto& a, const auto& w) {
                const auto ta = std::tan(t1);
                x -= x0;
                return a * ta * std::exp(-std::pow(x, 2) * std::pow(ta, 2) / w) *
                    std::exp(im * x * std::sin(t2) / std::sqrt(w));
            }, ca, cw
        );
    }

    template<typename Val, typename Arg, typename AV, typename WV>
    auto gaussian_source(const Arg& a, const Arg& b, const size_t& n, const Arg& x0, const AV& ca, const WV& cw) {
        return initial_conditions<Arg, Val>::create(
            [x0](const Arg& x, const auto& a, const auto& w) {
                return a * std::exp(-std::pow(x - x0, 2) / w);
            }, ca, cw
        );
    }

    template<typename Val, typename Arg, typename AV, typename WV>
    auto greene_source(const Arg& a, const Arg& b, const size_t& n, const Arg& x0, const AV& ca, const WV& cw) {
        return initial_conditions<Arg, Val>::create(
            [x0](const Arg& x, const auto& a, const auto& w) {
                const auto s = std::pow(x - x0, 2) / w;
                return a * (Arg(1.4467) - Arg(0.8402) * s) * std::exp(-s / 1.5256);
            }, ca, cw
        );
    }

    template<typename Arg, typename K0, typename PS, typename TA, typename Val = std::complex<Arg>>
    auto ray_source(const Arg& x,
                    const Arg& x0, const Arg& y0, 
                    const Arg& l1, const size_t& nl, 
                    const Arg& a0, const Arg& a1, const size_t& na,
                    const K0& k0,  const PS& ps, const utils::linear_interpolated_data_1d<Arg, Arg>& k_j,
                    const TA& tapering) {
        utils::dynamic_assert(k0.size() == k_j.size() && ps.size() == k_j.size(), "ray source: arguments k0, ps and k_j must have the same size");

        return initial_conditions<Arg, Val>(
            [=](const Arg& yl, const Arg& yr, const size_t& ny) {
                const auto [rx, ry] = rays::compute(x0, y0, l1, nl, a0, a1, na, k_j);
                const auto& as = rx.template get<0>();
                const auto& ls = rx.template get<1>();

                const auto ya = rays::_impl::calc_derivative_x(ry);

                size_t j = 0, i = 0;
                types::vector2d_t<size_t> li(k_j.size(), types::vector1d_t<size_t>(na));

                for (i = 0; i < as.size(); ++i) {
                    const auto k = std::min(size_t(l1 * std::cos(as[i]) / x), nl - 1);
                    for (j = 0; j < k_j.size(); ++j)
                        li[j][i] = x > rx[j][i][k]
                                   ? _impl::find_greater(x, rx[j][i], k)
                                   : _impl::find_less(x, rx[j][i], k);
                }

                const auto hl = l1 / (nl - 1);

                static constexpr auto ii = Val(0, 1);
                types::vector2d_t<std::tuple<Arg, Val>> pairs(k_j.size(), types::vector1d_t<std::tuple<Arg, Val>>(na));

                const auto lk = utils::function_as_vector([&i, &j, &k_j, &ry=ry](const size_t& l) { return k_j[j].point(ry[j][i][l]); }, 0);

                for (j = 0; j < k_j.size(); ++j) {
                    const auto m0 = std::exp(ii * Arg(M_PI) / Arg(4)) / std::sqrt(Arg(8) * Arg(M_PI) * k0[j]);

                    for (i = 0; i < na; ++i) {
                        const auto& ll = li[j][i];

                        const auto c = (x - rx[j][i][ll]) / (rx[j][i][ll + 1] - rx[j][i][ll]);
                        const auto ys = ry[j][i][ll] + c * (ry[j][i][ll + 1] - ry[j][i][ll]);

                        const auto lh = hl * c;
                        const auto sl = ls[ll] + lh;
                        const auto kl = k_j[j].point(ys);

                        const auto s = (utils::integrate_vector(ll + 1, hl, lk) + lh / Arg(2) * (kl + k_j[j].point(ry[j][i][ll]))) / k0[j];
                        const auto m = m0 / kl * k0[j] * std::sqrt(std::cos(as[i]) / ya[j].point(as[i], sl));

                        pairs[j][i] = std::make_tuple(ys, m * std::exp(ii * k0[j] * s));
                    }
                }

                const auto comp = [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); };
                for (j = 0; j < pairs.size(); ++j)
                    std::sort(pairs[j].begin(), pairs[j].end(), comp);

                const auto hy = (yr - yl) / (ny - 1);
                const auto yy = utils::mesh_1d(yl, yr, ny);

                types::vector1d_t<Arg> aa(ny, Arg(0));
                types::vector2d_t<Val> result(k_j.size(), types::vector1d_t<Val>(ny, Val(0)));

                const auto coords = utils::function_as_vector([&j, &pairs](const size_t& i) { return std::get<0>(pairs[j][i]); }, na);
                const auto values = utils::function_as_vector([&j, &pairs](const size_t& i) { return std::get<1>(pairs[j][i]); }, na);

                for (j = 0; j < k_j.size(); ++j) {
                    const auto il = std::min(size_t((std::get<0>(pairs[j][0]) - yl) / hy) + 1, ny - 1);
                    const auto ir = std::min(size_t((std::get<0>(pairs[j].back()) - yl) / hy), ny - 1);

                    utils::_impl::linear_interpolation::line(yy[il], yy[ir], coords, values, result[j].begin() + il, result[j].begin() + ir + 1);
                    utils::_impl::linear_interpolation::line(yy[il], yy[ir], coords, as, aa.begin() + il, aa.begin() + ir + 1);

                    tapering.interval(il, ir, yy, aa, result[j]);

                    std::transform(result[j].begin(), result[j].end(), result[j].begin(),
                                   [&ps, &j, k=std::exp(ii * k0[j] * x)](const auto& v) { return v * ps[j] / k; });
                }

                return result;
            }
        );
    }

    template<typename Arg, typename K0, typename PS, typename TA, typename Val = std::complex<Arg>>
    auto simple_ray_source(const Arg& x, const Arg& a0, const Arg& a1, const K0& k0,  const PS& ps, const TA& tapering) {
        utils::dynamic_assert(k0.size() == ps.size(), "simple ray source: arguments k0 and ps must have the same size");

        return initial_conditions<Arg, Val>(
            [=](const Arg& yl, const Arg& yr, const size_t& ny) {
                const auto ys = utils::mesh_1d(yl, yr, ny);
                const auto ls = utils::make_vector(ys, [x=std::pow(x, 2)](const auto& y) { return std::sqrt(y * y + x); });
                const auto sl = utils::make_vector(ls, [](const auto& l) { return std::sqrt(l); });
                const auto as = utils::make_vector(ys, [&x](const auto& y) { return std::atan2(y, x); });

                static constexpr auto ii = Val(0, 1);

                types::vector2d_t<Val> result(k0.size(), types::vector1d_t<Val>(ny));
                for (size_t j = 0; j < k0.size(); ++j) {
                    const auto m0 = std::exp(ii * (Arg(M_PI) / Arg(4) - k0[j] * x)) / std::sqrt(Arg(8) * Arg(M_PI) * k0[j]);
                    for (size_t i = 0; i < ny; ++i)
                        result[j][i] = ps[j] * m0 * std::exp(ii * k0[j] * ls[i]) / sl[i];
                }

                const auto il = std::distance(as.begin(), std::upper_bound(as.begin(), as.end(), a0));
                const auto ir = std::distance(as.begin(), std::lower_bound(as.begin(), as.end(), a1));

                tapering.tail_all(il, ir, ys, as, result);

                return result;
            }
        );
    }

}// namespace ample
