#pragma once
#include <tuple>
#include <vector>
#include <complex>
#include <istream>
#include <exception>
#include <type_traits>

namespace acstc {

    namespace types {

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

        using rvector1d_t = vector1d_t<real_t>;
        using rvector2d_t = vector2d_t<real_t>;
        using cvector1d_t = vector1d_t<complex_t>;
        using cvector2d_t = vector2d_t<complex_t>;

        template<typename T>
        using tuple2_t = std::tuple<T, T>;

        template<typename T>
        struct vector_dim {

            static constexpr size_t value = 0;            

        };

        template<typename T>
        struct vector_dim<std::vector<T>> {

            static constexpr size_t value = vector_dim<T>::value + 1;

        };

        template<typename T>
        constexpr size_t vector_dim_v = vector_dim<std::decay_t<T>>::value;

    }// namespace types

}// namespace acstc
