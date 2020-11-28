#pragma once
#include <string>
#include <stdexcept>

namespace acstc::utils {

    template<typename E = std::runtime_error>
    void dynamic_assert(const bool& condition, const char* message) {
        if (!condition)
            throw E(message);
    }

    template<typename E = std::runtime_error>
    void dynamic_assert(const bool& condition, const std::string& message) {
        dynamic_assert(condition, message.c_str());
    }

    template<typename E = std::runtime_error>
    void dynamic_assert(const bool& condition) {
        dynamic_assert<E>(condition, "Dynamic assertion failed");
    }

}// namespace acstc
