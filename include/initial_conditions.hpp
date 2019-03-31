#pragma once

#include <cmath>
#include <tuple>
#include <cstddef>
#include <functional>
#include "utils/types.hpp"
#include "utils/utils.hpp"

namespace acstc {

    template<typename Arg, typename Val = Arg>
    class initial_conditions {

    public:

        initial_conditions() = delete;
        initial_conditions(const initial_conditions&) = delete;
        initial_conditions(initial_conditions&&) = delete;
        ~initial_conditions() = delete;

        template<typename F, typename... Args>
        static auto create(const Arg& a, const Arg& b, const size_t& n, const F& func, const Args&... args) {
            const auto h = (b - a) / (n - 1);
            types::vector1d_t<Arg> xs(n);
            for (size_t i = 0; i < n; ++i)
                xs[i] = a + i * h;
            return create(xs, func, std::forward<const Args>(args)...);
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

}// namespace acstc
