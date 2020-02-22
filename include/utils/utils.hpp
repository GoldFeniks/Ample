#pragma once
#include <tuple>
#include <string>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <filesystem>
#include "types.hpp"

namespace acstc {

    namespace utils {

        namespace __impl {

            template<size_t N, size_t M = N>
            struct gather_values {

                gather_values() = delete;

                template<typename Tuple, typename... Values>
                static auto get(const size_t& index, const Tuple& tuple, Values&... values) {
                    return gather_values<N, M - 1>::get(
                            index, tuple, std::forward<Values>(values)..., std::get<N - M>(tuple)[index]);
                }

            };

            template<size_t N>
            struct gather_values<N, 0> {

                gather_values() = delete;

                template<typename Tuple, typename... Values>
                static auto get(const size_t& index, const Tuple& tuple, Values&... values) {
                    return std::make_tuple(std::forward<const Values>(values)...);
                }

            };

            template<typename It>
            class stride_iterator {

            public:

                explicit stride_iterator(It begin, It end, const size_t& k) : _it(begin), _end(end), _k(k) {}

                bool operator==(const stride_iterator& other) const {
                    return _it == other._it;
                }

                bool operator!=(const stride_iterator& other) const {
                    return !(*this == other);
                }

                bool operator==(const It& other) const {
                    return _it == other;
                }

                bool operator!=(const It& other) const {
                    return !(*this == other);
                }

                const auto& operator*() const {
                    return *_it;
                }

                stride_iterator& operator++() {
                    advance();
                    return *this;
                }

                const stride_iterator operator++(int) {
                    auto res = stride_iterator(_it, _end, _k);
                    advance();
                    return res;
                }

             private:

                void advance() {
                    if (std::distance(_it, _end) < _k)
                        _it = _end;
                    else
                        _it += _k;
                }

                It _it, _end;
                const size_t _k;

            };

            template<size_t>
            struct __integrator;

            template<>
            struct __integrator<5> {

                template<typename T, typename F>
                static T integrate(const T& x, const T& h, const F& func, T& fx) {
                    return T(2) * h / T(45) * (
                        T(7)  * (fx + (fx = func(x + T(4) * h))) + 
                        T(32) * (func(x + h) + func(x + T(3) * h)) + 
                        T(12) * func(x + T(2) * h)
                    );
                }

                template<typename T, typename V>
                static T integrate_vector(const size_t& i, const T& h, const V& data, T& fx) {
                    return T(2) * h / T(45) * (
                        T(7)  * (fx + (fx = data[i + 4])) +
                        T(32) * (data[i + 1] + data[i + 3]) +
                        T(12) * data[i + 2]
                    );
                }

            };

            template<>
            struct __integrator<4> {

                template<typename T, typename F>
                static T integrate(const T& x, const T& h, const F& func, T& fx) {
                    return T(3) * h / T(8) * (fx + T(3) * (func(x + h) + func(x * T(2) * h)) + (fx = func(x + T(3) * h)));
                }

                template<typename T, typename V>
                static T integrate_vector(const size_t& i, const T& h, const V& data, T& fx) {
                    return T(3) * h / T(8) * (fx + T(3) * (data[i + 1] + data[i + 2]) + (fx = data[i + 3]));
                }

            };

            template<>
            struct __integrator<3> {

                template<typename T, typename F>
                static T integrate(const T& x, const T& h, const F& func, T& fx) {
                    return h / T(3) * (fx + T(4) * func(x + h) + (fx = func(x + T(2) * h)));
                }

                template<typename T, typename V>
                static T integrate_vector(const size_t& i, const T& h, const V& data, T& fx) {
                    return h / T(3) * (fx + T(4) * data[i + 1] + (fx = data[i + 2]));
                }

            };

            template<>
            struct __integrator<2> {

                template<typename T, typename F>
                static T integrate(const T& x, const T& h, const F& func, T& fx) {
                    return h / T(2) * (fx + (fx = func(x + h)));
                }

                template<typename T, typename V>
                static T integrate_vector(const size_t& i, const T& h, const V& data, T& fx) {
                    return h / T(2) * (fx + (fx = data[i + 1]));
                }

            };

            template<size_t N>
            struct integrator {

                template<typename T, typename F>
                static T integrate(T x, const T& h, const size_t& n, const F& func) {
                    T res = T(0), fx = func(x);
                    for (size_t i = 0; i < n; ++i, x += (N - 1) * h)
                        res += __integrator<N>::integrate(x, h, func, fx);
                    return res;
                }

                template<typename T, typename F>
                static T integrate_all(const T& a, const T& h, const size_t& n, const F& func) {
                    const auto k = n / (N - 1);
                    const auto res = k ? integrate(a, h, k, func) : T(0);
                    if constexpr (N > 2)
                        if (const auto m = n - k * (N - 1))
                            return res + integrator<N - 1>::integrate_all(a + k * (N - 1) * h, h, m, func);
                    return res;
                }

                template<typename T, typename V>
                static T integrate_vector(size_t a, const T& h, const size_t& n, const V& data) {
                    T res = T(0), fx = data[a];
                    for (size_t i = 0; i < n; ++i, a += (N - 1))
                        res += __integrator<N>::integrate_vector(a, h, data, fx);
                    return res;
                }

                template<typename T, typename V>
                static T integrate_vector_all(const size_t& a, const T& h, const size_t& n, const V& data) {
                    const auto k = n / (N - 1);
                    const auto res = k ? integrate_vector(a, h, k, data) : T(0);
                    if constexpr (N > 2)
                        if (const auto m = n - k * (N - 1))
                            return res + integrator<N - 1>::integrate_vector_all(a + k * (N - 1), h, m, data);
                    return res;
                }

            };

        }// namespace __impl

        template<typename... T>
        class zip {

        public:

            explicit zip(const T&... values)
                : _values(std::forward<const T>(values)...), _size(std::get<0>(_values).size()) {}

            auto begin() const {
                return iterator(*this);
            }

            auto end() const {
                return iterator(*this, size());
            }

            auto size() const {
                return _size;
            }

            auto operator[](const size_t& index) const {
                return __impl::gather_values<N>::get(index, _values);
            }

        private:

            const std::tuple<const T&...> _values;
            const size_t _size;

            static constexpr auto N = sizeof...(T);

            class iterator {

            public:

                explicit iterator(const zip& owner, const size_t index = 0) : _owner(owner), _index(index) {}

                bool operator==(const iterator& other) const {
                    return _index == other._index;
                }

                bool operator!=(const iterator& other) const {
                    return !(*this == other);
                }

                auto operator*() const {
                    return _owner[_index];
                }

                iterator& operator++() {
                    ++_index;
                    return *this;
                }

                iterator operator++(int) {
                    auto res = iterator(_owner, _index);
                    ++_index;
                    return res;
                }

            private:

                const zip& _owner;
                size_t _index = 0;

            };

        };

        template<typename F>
        class function_as_vector {

        public:

            function_as_vector(const F& func, const size_t& size) : _func(func), _size(size) {}

            const auto operator[](const size_t& i) const {
                return _func(i);
            }

            const size_t& size() const {
                return _size;
            }

        private:

            const F _func;
            const size_t _size;

        };

        template<typename V, typename T>
        auto find_indices(const V& values, const T& value) {
            auto it = std::lower_bound(values.begin(), values.end(), value);
            if (it == values.end())
                return std::tuple<size_t, size_t>(0, 1);
            if (it == values.begin())
                ++it;
            const auto d = std::distance(values.begin(), it);
            return std::tuple<size_t, size_t>(d - 1, d);
        }

        template<typename T>
        auto mesh_1d(const T& a, const T& b, const size_t& n) {
            const auto h = (b - a) / (n - 1);
            types::vector1d_t<T> result;
            result.reserve(n);
            for (size_t i = 0; i < n; ++i)
                result.emplace_back(a + i * h);
            return result;
        }

        std::filesystem::path make_file_path(const std::filesystem::path root, const std::filesystem::path& filename) {
            if (!root.has_filename() || filename.is_absolute())
                return filename;
            std::filesystem::path result = root;
            return result.replace_filename(filename);
        }

        template<typename It>
        auto stride(It begin, It end, const size_t& k) {
            return std::make_tuple(
                __impl::stride_iterator(begin, end, k),
                __impl::stride_iterator(end,   end, k)
            );
        }

        template<typename T, typename F>
        T integrate(const T& a, const T& b, const size_t& n, const F& func) {
            const auto& h = (b - a) / (n - 1);
            return __impl::integrator<5>::integrate_all(a, h, n, func);
        }

        template<typename T, typename V>
        T integrate_vector(const size_t& a, const size_t& b, const T& h, const V& data) {
            return __impl::integrator<5>::integrate_vector_all(a, h, b - a - 1, data);
        }

        template<typename T, typename V>
        T integrate_vector(const size_t& b, const T& h, const V& data) {
            return integrate_vector(size_t(0), b, h, data);
        }

    }// namespace utils

}// namespace acstc