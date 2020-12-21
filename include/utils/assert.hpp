#pragma once
#include <stdexcept>
#include "join.hpp"

namespace acstc::utils {

    template<typename E>
    struct set_exception {

        template<typename... T>
        static void dynamic_assert(const bool& condition, const T&... args) {
            if (!condition)
                throw E(join(args...));
        }

    };

    template<typename... T>
    void dynamic_assert(const bool& condition, const T&... args) {
        set_exception<std::runtime_error>::dynamic_assert(condition, args...);
    }

    template<>
    void dynamic_assert<>(const bool& condition) {
        dynamic_assert(condition, "Dynamic assertion failed");
    }

}// namespace acstc
