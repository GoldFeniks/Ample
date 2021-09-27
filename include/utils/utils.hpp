#pragma once
#include <array>
#include <tuple>
#include <string>
#include <sstream>
#include <cstddef>
#include <iterator>
#include <optional>
#include <algorithm>
#include <filesystem>
#include <functional>
#include <type_traits>
#include "types.hpp"

namespace ample::utils {

    namespace _impl {

        template<typename It>
        class stride_iterator {

        public:

            using iterator_category = std::output_iterator_tag;
            using value_type        = typename It::value_type;
            using difference_type   = typename It::difference_type;
            using pointer           = typename It::pointer;
            using reference         = typename It::reference;

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

            decltype(*std::declval<It>()) operator*() const {
                return *_it;
            }

            decltype(*std::declval<It>()) operator->() const {
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
        struct _integrator;

        template<>
        struct _integrator<5> {

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
        struct _integrator<4> {

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
        struct _integrator<3> {

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
        struct _integrator<2> {

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
                    res += _integrator<N>::integrate(x, h, func, fx);
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
                    res += _integrator<N>::integrate_vector(a, h, data, fx);
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

    }// namespace _impl

    template<typename F>
    class function_as_vector {

    public:

        function_as_vector(const F& func, const size_t& size) : _func(func), _size(size) {}

        auto operator[](const size_t& i) const {
            return _func(i);
        }

        [[nodiscard]] const size_t& size() const {
            return _size;
        }

        auto front() const {
            return _func(0);
        }

        auto back() const {
            return _func(_size - 1);
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
    auto mesh_1d(const T& a, const T& b, const size_t& n) -> types::vector1d_t<decltype(b - a)> {
        if (n == 1)
            return types::vector1d_t<T>{ a };

        const auto h = (b - a) / static_cast<T>(n - 1);
        types::vector1d_t<T> result;
        result.reserve(n);
        for (size_t i = 0; i < n; ++i)
            result.emplace_back(a + static_cast<T>(i) * h);
        return result;
    }

    template<typename T>
    void mesh_1d(...);

    template<typename V, typename F, typename = std::enable_if_t<std::is_invocable_v<F, decltype(std::declval<V>()[0])>>>
    auto make_vector(const V& data, F&& func) {
        types::vector1d_t<decltype(func(data[0]))> result(data.size());
        std::transform(data.begin(), data.end(), result.begin(), func);
        return result;
    }

    template<typename V, typename F, typename = std::enable_if_t<std::is_invocable_v<F, decltype(std::declval<V>()[0]), size_t>>>
    auto make_vector_i(const V& data, F&& func) {
        return make_vector(data, [&func, i=size_t(0)](const auto& value) mutable { return func(value, i++); });
    }

    template<typename F, typename = std::enable_if_t<std::is_invocable_v<F, size_t>>>
    auto make_vector_i(const size_t& n, F&& func) {
        types::vector1d_t<decltype(func(size_t(0)))> result;
        for (size_t i = 0; i < n; ++i)
            result.emplace_back(func(i));
        return result;
    }

    std::filesystem::path make_file_path(const std::filesystem::path& root, const std::filesystem::path& filename) {
        if (!root.has_filename() || filename.is_absolute())
            return filename;
        std::filesystem::path result = root;
        return result.replace_filename(filename);
    }

    template<typename It>
    auto stride(It begin, It end, const size_t& k) {
        return std::make_tuple(
                _impl::stride_iterator(begin, end, k),
                _impl::stride_iterator(end, end, k)
        );
    }

    template<typename T, typename F>
    T integrate(const T& a, const T& b, const size_t& n, const F& func) {
        const auto& h = (b - a) / (n - 1);
        return _impl::integrator<5>::integrate_all(a, h, n, func);
    }

    template<typename T, typename V>
    T integrate_vector(const size_t& a, const size_t& b, const T& h, const V& data) {
        return _impl::integrator<5>::integrate_vector_all(a, h, b - a - 1, data);
    }

    template<typename T, typename V>
    T integrate_vector(const size_t& b, const T& h, const V& data) {
        return integrate_vector(size_t(0), b, h, data);
    }

    template<typename T>
    class span {

    public:

        span() = default;
        span(T* values, const size_t& n) : _n(n), _values(values), _has_value(true) {}

        void assign(T* values, const size_t& n) {
            _n = n;
            _values = values;
            _has_value = true;
        }

        auto& operator[](const size_t& index) {
            return _values[index];
        }

        const auto& operator[](const size_t& index) const {
            return _values[index];
        }

        [[nodiscard]] size_t size() const {
            return _n;
        }

        T* data() {
            return _values;
        }

        const T* data() const {
            return _values;
        }

        [[nodiscard]] bool has_value() const {
            return _has_value;
        }

        void reset() {
            _n = 0;
            _values = nullptr;
            _has_value = false;
        }

    private:

        size_t _n = 0;
        T* _values = nullptr;
        bool _has_value = false;

    };

    template<typename T>
    class lazy_value {

    public:

        lazy_value(std::function<T()> getter) : _getter(std::move(getter)) {}

        const T& get() const {
            if (!_value.has_value())
                _value = _getter();
            return _value.value();
        }

        const T& operator()() const {
            return get();
        }

    private:

        mutable std::optional<T> _value;
        std::function<T()> _getter;

    };

}// namespace ample::utils
