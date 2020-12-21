#pragma once
#include <tuple>
#include <vector>
#include <complex>
#include <istream>
#include <exception>
#include <type_traits>

namespace acstc::types {

    template<typename T, size_t N>
    struct vectornd {

        using type = std::vector<typename vectornd<T, N - 1>::type>;

    };

    template<typename T>
    struct vectornd<T, 0> {

        using type = T;

    };

    template<typename T, size_t N>
    using vectornd_t = typename vectornd<T, N>::type;

    using real_t = double;
    using complex_t = std::complex<real_t>;

    template<typename T>
    using vector1d_t = std::vector<T>;

    template<typename T>
    using vector2d_t = std::vector<std::vector<T>>;

    template<typename T>
    using vector3d_t = std::vector<std::vector<std::vector<T>>>;

    template<typename T>
    using vector4d_t = std::vector<std::vector<std::vector<std::vector<T>>>>;

    template<typename T>
    using tuple2_t = std::tuple<T, T>;

    template<typename T>
    struct point { T x, y, z; };

}// namespace acstc::types
