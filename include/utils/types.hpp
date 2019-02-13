#pragma once
#include <vector>
#include <complex>

namespace acstc {

    namespace types {

        using real_t = double;
        using complex_t = std::complex<real_t>;

        template<typename T>
        using vector1d_t = std::vector<T>;

        template<typename T>
        using vector2d_t = std::vector<std::vector<T>>;

        using cvector1d_t = vector1d_t<complex_t>;
        using cvector2d_t = vector2d_t<complex_t>;

    }// namespace types

}// namespace acstc
