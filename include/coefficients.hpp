#pragma once
#include <cmath>
#include <tuple>
#include <cstddef>
#include <utility>
#include <algorithm>
#include <functional>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/assert.hpp"
#include "utils/comparators.hpp"

namespace ample {

    namespace _impl {

        template<typename T>
        using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;

        template<typename T>
        using matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

        template<typename T>
        struct solver_type {

            using type = Eigen::EigenSolver<matrix_t<T>>;

        };

        template<typename T>
        struct solver_type<std::complex<T>> {

            using type = Eigen::ComplexEigenSolver<matrix_t<std::complex<T>>>;

        };


        template<typename T>
        using solver_type_t = typename solver_type<T>::type;

        template<typename T>
        struct pade_coefficients {

            auto operator()(const types::vector1d_t<T>& c, const size_t& n, const size_t& m) {
                vector_t<T> b(m);
                matrix_t<T> A(m, m);

                for (size_t i = 0, ci = n + 1; i < m; ++i, ++ci) {
                    b(i) = -c[ci];

                    for (size_t j = 0; j < std::min(ci, m); ++j)
                        A(i, j) = c[ci - j - 1];
                }

                const vector_t<T> x = A.colPivHouseholderQr().solve(b);

                types::vector1d_t<T> ac(n + 1), bc(m + 1, T(1));
                for (size_t i = 0; i < n + 1; ++i) {
                    ac[i] = c[i];
                    for (size_t j = 0; j < std::min(i, m); ++j)
                        ac[i] += c[i - j - 1] * x(j);
                }

                for (size_t i = 0; i < m; ++i)
                    bc[i + 1] = x(i);

                return std::make_tuple(ac, bc);
            }

        };

        template<typename T>
        struct theta_pade_coefficients {

            T theta{};

            auto operator()(const types::vector1d_t<T>& c, const size_t& n, const size_t& m) {
                utils::dynamic_assert(n == m, "coefficients: n(", n, ") must be equal to m(", m, ") for theta-weighted coefficients");

                pade_coefficients<T> coefficients{};
                const auto [a0, b0] = coefficients(c, n - 1, n);
                      auto [a1, b1] = coefficients(c, n, n);

                const auto theta1 = T(1) - theta;
                for (size_t i = 0; i < a0.size(); ++i)
                    a1[i] = a1[i] * theta + a0[i] * theta1;

                a1.back() *= theta;

                for (size_t i = 0; i < b0.size(); ++i)
                    b1[i] = b1[i] * theta + b0[i] * theta1;

                return std::make_tuple(a1, b1);
            }

        };

        template<typename T>
        struct exp_taylor {

            static constexpr auto on = T(1);
            static constexpr auto tw = T(2);
            static constexpr auto fo = T(4);
            static constexpr auto si = T(6);

            auto operator()(const T& value, const size_t& n) {
                const auto v2 = std::pow(value, 2);

                types::vector1d_t<T> c(n);
                c[0] = on;
                c[1] = value /tw;

                for (size_t i = 2, m = 0; i < n; ++i, ++m) {
                    const auto mm = T(m);
                    c[i] = (v2 * c[m] - (tw + si * mm + fo * mm * mm) * c[m + 1]) / (fo * (on + mm) * (tw + mm));
                }

                return c;
            }

        };

        template<typename T>
        struct root_taylor {

            auto operator()(const T& value, const size_t& n) {
                types::vector1d_t<T> c(n, T(0));
                std::for_each(c.begin() + 1, c.end(), 
                    [&value, p=T(1), i=size_t(0), f=T(1)](auto& x) mutable {
                        ++i;
                        x = value * (p *= (T(1.5) - T(i))) / (f *= i);
                    }
                );

                return c;
            }

        };

        template<typename T>
        struct no_transform {

            auto operator()(types::vector1d_t<T>& p, types::vector1d_t<T>& q) {}

        };

        template<typename T>
        struct wampe_transform {

            auto operator()(types::vector1d_t<T>& p, types::vector1d_t<T>& q) {
                for (size_t i = 0; i < p.size(); ++i) {
                    const auto q2 = T(2) * q[i];
                    std::tie(p[i], q[i]) = std::make_tuple(-p[i] - q2, p[i] - q2);
                }

                for (size_t i = p.size(); i < q.size(); ++i) {
                    const auto q2 = T(-2) * q[i];
                    p.push_back(q2);
                    q[i] = q2;
                }
            }

        };

        template<typename T>
        auto get_roots(const types::vector1d_t<T>& c) {
            const auto n = c.size() - 1;
            matrix_t<T> A = matrix_t<T>::Zero(n, n);

            for (size_t i = 0; i < n - 1; ++i)
                A(i + 1, i) = T(1);

            for (size_t i = 0; i < n; ++i)
                A(0, i) = -c[n - i - 1] / c.back();

            solver_type_t<T> solver;
            solver.compute(A);
            const vector_t<T> values = solver.eigenvalues();

            auto result = utils::make_vector_i(n, [&values](const auto& i) { return values(i); });
            std::sort(result.begin(), result.end(), utils::less<T>{});
            return result;
        }

    }// namespace ample::_impl

    template<typename T>
    class coefficients {

    public:

        using taylor_func_t      = std::function<types::vector1d_t<T>(const T&, const size_t&)>;
        using transform_func_t   = std::function<void (types::vector1d_t<T>&, types::vector1d_t<T>&)>;
        using coeficients_func_t = std::function<std::tuple<types::vector1d_t<T>, types::vector1d_t<T>>(const types::vector1d_t<T>&, const size_t&, const size_t&)>;

        coefficients(const size_t& n,
                     taylor_func_t taylor,
                     transform_func_t transform,
                     coeficients_func_t coefficients) : ample::coefficients<T>(n, n, std::move(taylor), std::move(transform), std::move(coefficients)) {}

        coefficients(const size_t& n, const size_t& m,
                     taylor_func_t taylor,
                     transform_func_t transform,
                     coeficients_func_t coefficients) : _n(n), _m(m), _taylor(std::move(taylor)), _transform(std::move(transform)), _coefficients(std::move(coefficients)) {
            utils::dynamic_assert(n > 0,  "coefficients: n(", n, ") must be positive");
            utils::dynamic_assert(n <= m, "coefficients: n(", n, ") must be less or equal to m(", m, ")");
        }

        auto get(const T& value) const {
            const auto tc = _taylor(value, _n + _m + 1);
            auto [np, dp] = _coefficients(tc, _n, _m);

            _transform(np, dp);

            const auto nr = _impl::get_roots<T>(np);
            const auto dr = _impl::get_roots<T>(dp);

            types::vector1d_t<size_t> ix(_m);
            std::iota(ix.begin(), ix.end(), size_t(0));

            auto a = types::vector1d_t<T>(_m);
            for (size_t i = 0; i < _m; ++i) {
                const auto p = std::accumulate(nr.begin(), nr.end(), np.back(), [&x=dr[i]](const auto& a, const auto& x0) { return a * (x - x0); });
                const auto q = std::accumulate(ix.begin(), ix.end(), dp.back(),
                                               [&x=dr[i], &i, &dr](const auto& a, const auto& j) {
                                                   return i == j ? a : a * (x - dr[j]);
                                               }
                );

                a[i] = -p / q / dr[i];
            }

            auto a0 = np.front() / dp.front() - std::accumulate(a.begin(), a.end(), T(0));
            return std::make_tuple(a0, std::move(a), utils::make_vector(dr, [](const auto& x) { return -T(1) / x; }));
        }

        [[nodiscard]] auto nc() const {
            return _m;
        }

    private:

        size_t _n{}, _m{};
        taylor_func_t _taylor;
        transform_func_t _transform;
        coeficients_func_t _coefficients;

    };

    template<typename T>
    auto ssp_coefficients(const size_t& n, const size_t& m) {
        return coefficients<T>(n, m, _impl::exp_taylor<T>{}, _impl::no_transform<T>{}, _impl::pade_coefficients<T>{});
    }

    template<typename T>
    auto ssp_coefficients(const size_t& n) {
        return ssp_coefficients<T>(n, n);
    }

    template<typename T>
    auto theta_ssp_coefficients(const size_t& n, const T& theta) {
        return coefficients<T>(n, n, _impl::exp_taylor<T>{}, _impl::no_transform<T>{}, _impl::theta_pade_coefficients<T>{ theta });
    }

    template<typename T>
    auto wampe_coefficients(const size_t& n, const size_t& m) {
        return coefficients<T>(n, m, _impl::root_taylor<T>{}, _impl::wampe_transform<T>{}, _impl::pade_coefficients<T>{});
    }

    template<typename T>
    auto wampe_coefficients(const size_t& n) {
        return wampe_coefficients<T>(n, n);
    }

    template<typename T>
    auto theta_wampe_coefficients(const size_t& n, const T& theta) {
        return coefficients<T>(n, n, _impl::root_taylor<T>{}, _impl::wampe_transform<T>{}, _impl::theta_pade_coefficients<T>{ theta });
    }

}// namespace ample
