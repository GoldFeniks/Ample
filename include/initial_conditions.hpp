#pragma once

#include <cmath>
#include <tuple>
#include <complex>
#include <cstddef>
#include <algorithm>
#include <functional>
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/interpolation.hpp"

namespace acstc {

    namespace __impl {

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
        inline auto smooth(const size_t& l, const size_t& r, const T& rc, const T& d, const C& coords, V& result) {
            for (size_t i = l + 1; i < r; ++i) {
                const auto v = std::pow(rc - coords[i], 2) / d;
                result[i] *= std::exp(v * v * (T(1) / (v - T(1)) - T(1)));
            }
        }

    }// namespace __impl
    

    template<typename Arg, typename Val = Arg>
    class initial_conditions {

    public:

        initial_conditions() = delete;
        initial_conditions(const initial_conditions&) = delete;
        initial_conditions(initial_conditions&&) = delete;
        ~initial_conditions() = delete;

        template<typename F, typename... Args>
        static auto create(const Arg& a, const Arg& b, const size_t& n, const F& func, const Args&... args) {
            return create(utils::mesh_1d(a, b, n), func, std::forward<const Args>(args)...);
        }

        template<typename V, typename F, typename... Args>
        static auto create(const V& xs, const F& func, const Args&... args) {
            const auto zip = utils::zip(std::forward<const Args>(args)...);
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
        return initial_conditions<Arg, Val>::create(a, b, n,
                [x0, t1, t2, im=Val(0, 1)](Arg x, const auto& a, const auto& w) {
                    const auto ta = std::tan(t1);
                    x -= x0;
                    return a * ta * std::exp(-std::pow(x, 2) * std::pow(ta, 2) / w) *
                        std::exp(im * x * std::sin(t2) / std::sqrt(w));
                }, ca, cw);
    }

    template<typename Val, typename Arg, typename AV, typename WV>
    auto gaussian_source(const Arg& a, const Arg& b, const size_t& n, const Arg& x0, const AV& ca, const WV& cw) {
        return initial_conditions<Arg, Val>::create(a, b, n,
                [x0](const Arg& x, const auto& a, const auto& w) {
                    return a * std::exp(-std::pow(x - x0, 2) / w);
                }, ca, cw);
    }

    template<typename Val, typename Arg, typename AV, typename WV>
    auto greene_source(const Arg& a, const Arg& b, const size_t& n, const Arg& x0, const AV& ca, const WV& cw) {
        return initial_conditions<Arg, Val>::create(a, b, n,
                [x0](const Arg& x, const auto& a, const auto& w) {
                    const auto s = std::pow(x - x0, 2) / w;
                    return a * (Arg(1.4467) - Arg(0.8402) * s) * std::exp(-s / 1.5256);
                }, ca, cw);
    }

    template<typename Arg, typename K0, typename PS, typename Val = std::complex<Arg>>
    auto ray_source(const Arg& x, 
                    const Arg& yl, const Arg& yr, const size_t& ny,
                    const Arg& x0, const Arg& y0, 
                    const Arg& l1, const size_t& nl, 
                    const Arg& a0, const Arg& a1, const size_t& na,
                    const K0& k0,  const PS& ps, const utils::linear_interpolated_data_1d<Arg, Arg>& k_j,
                    const double smooth_ratio = 0.1) {
        if (k0.size() != k_j.size() || ps.size() != k_j.size())
            throw std::logic_error("Arguments k0, ps and k_j must have the same size");

        const auto [rx, ry] = rays::compute(x0, y0, l1, nl, a0, a1, na, k_j);
        const auto& as = rx.template get<0>();
        const auto& ls = rx.template get<1>();

        const auto ya = rays::__impl::calc_derivative_x(ry);

        size_t j = 0, i = 0;
        types::vector2d_t<size_t> li(k_j.size(), types::vector1d_t<size_t>(na));

        for (i = 0; i < as.size(); ++i) {
            const auto k = std::min(size_t(l1 * std::cos(as[i]) / x), nl - 1);
            for (j = 0; j < k_j.size(); ++j)
                li[j][i] = x > rx[j][i][k]
                    ? __impl::find_greater(x, rx[j][i], k)
                    : __impl::find_less(x, rx[j][i], k);
        }

        const auto hl = l1 / (nl - 1);

        static constexpr auto ii = Val(0, 1);
        types::vector2d_t<std::tuple<Arg, Val>> pairs(k_j.size(), types::vector1d_t<std::tuple<Arg, Val>>(na));

        const auto lk = utils::function_as_vector([&i, &j, &k_j, &ry](const size_t& l) { return k_j[j].point(ry[j][i][l]); }, 0);

        for (j = 0; j < k_j.size(); ++j) {
            const auto m0 = std::exp(ii * M_PI / Arg(4)) / std::sqrt(Arg(8) * M_PI * k0[j]);

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
        types::vector2d_t<Val> result(k_j.size(), types::vector1d_t<Val>(ny, Val(0)));        

        const auto coords = utils::function_as_vector([&j, &pairs](const size_t i) { return std::get<0>(pairs[j][i]); }, na);
        const auto values = utils::function_as_vector([&j, &pairs](const size_t i) { return std::get<1>(pairs[j][i]); }, na);

        for (j = 0; j < result.size(); ++j) {
            const auto il = std::min(size_t((std::get<0>(pairs[j][0]) - yl) / hy) + 1, ny - 1);
            const auto ir = std::min(size_t((std::get<0>(pairs[j].back()) - yl) / hy), ny - 1);
            const auto id = ir - il;
            const auto jl = il + size_t(id * smooth_ratio);
            const auto jr = ir - size_t(id * smooth_ratio);

            const auto ly = yy[jl];
            const auto ry = yy[jr];

            utils::__impl::linear_interpolation::line(ly, ry, coords, values, result[j].begin() + il, result[j].begin() + ir + 1);

            __impl::smooth(il, jl, ly, std::pow(ly - yy[il], 2), yy, result[j]);
            __impl::smooth(jr, ir, ry, std::pow(ry - yy[ir], 2), yy, result[j]);

            for (i = il + 1; i < ir; ++i)
                result[j][i] *= ps[j];

            result[j][il] = result[j][ir] = Val(0);
        }

        return result;        
    }

}// namespace acstc
