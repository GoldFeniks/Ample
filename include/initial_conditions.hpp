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

        template<typename F, typename... Args>
        explicit initial_conditions(const Arg& a, const Arg& b, const size_t& n, const F& func, const Args&... args) {
            const auto h = (b - a) / (n - 1);
            types::vector1d_t<Arg> xs(n);
            for (size_t i = 0; i < n; ++i)
                xs[i] = a + i * h;
            _init(func, xs, std::forward<const Args>(args)...);
        }

        template<typename V, typename F, typename... Args>
        explicit initial_conditions(const V& xs, const F& func, const Args&... args) {
            _init(func, xs, std::forward<const Args>(args)...);
        }

        auto size() const {
            return _data.size();
        }

        auto operator[](const size_t& index) const {
            return _data[index];
        }

        auto begin() const {
            return _data.begin();
        }

        auto end() const {
            return _data.end();
        }

    private:

        types::vector2d_t<Val> _data;

        template<typename F, typename... Args>
        void _init(const F& func, const types::vector1d_t<Arg>& xs, const Args&... args) {
            const auto zip = utils::zip(std::forward<const Args>(args)...);
            _data.resize(zip.size());
            for (size_t i = 0; i < size(); ++i) {
                _data[i].resize(xs.size());
                auto tuple = std::tuple_cat(std::make_tuple(xs[0]), zip[i]);
                for (size_t j = 0; j < xs.size(); ++j) {
                    std::get<0>(tuple) = xs[j];
                    _data[i][j] = std::apply(func, tuple);
                }
            }
        }

    };

    template<typename Arg, typename Val = Arg, typename KV>
    auto generalized_gaussian_source(const Arg& a, const Arg& b, const size_t& n, const Arg& a0,
            const Arg& t1, const Arg& t2, const KV& k) {
        return initial_conditions<Arg, Val>(a, b, n,
            [a0, t1, t2](Arg a, const auto& k) {
                const auto ta = std::tan(t1);
                a -= a0;
                return std::sqrt(k) * ta *
                    std::exp(-std::pow(k, 2) * std::pow(a, 2) * std::pow(ta, 2) / Arg(2)) *
                    std::exp(Val(0, 1) * k * a * std::sin(t2));
            }, k);
    }

    template<typename Arg, typename Val = Arg, typename KV>
    auto gaussian_source(const Arg& a, const Arg& b, const size_t& n, const Arg& a0, const KV& k) {
        return generalized_gaussian_source<Arg, Val, KV>(a, b, n, a0, Arg(M_PI / 4), Arg(0), k);
    }

    template<typename Arg, typename Val = Arg, typename KV>
    auto greene_source(const Arg& a, const Arg& b, const size_t& n, const Arg& a0, const KV& k) {
        return initial_conditions<Arg, Val>(a, b, n,
            [a0](const Arg& a, const auto& k) {
                const auto& s = std::pow(a - a0, 2) * std::pow(k, 2);
                return std::sqrt(k) * (Arg(1.4467) - Arg(0.4201) * s) * std::exp(-s / Arg(3.0512));
            }, k);
    }

}// namespace acstc
