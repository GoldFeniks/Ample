#pragma once
#include "types.hpp"

namespace ample::utils {

    template<typename T>
    struct less {

        bool operator()(const types::real_t& a, const types::real_t& b) {
            return a < b;
        }

    };

    template<typename T>
    struct less<std::complex<T>> {

        bool operator()(const std::complex<T>& a, const std::complex<T>& b) {
            return a.real() == b.real() ? a.imag() < b.imag() : a.real() < b.real();
        }

    };

}// namespace ample::utils
