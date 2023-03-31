#pragma once
#include <ranges>
#include <type_traits>

namespace ample::utils::concepts {

    template<typename V, typename T>
    concept vector1d = requires(V v, size_t n) {
        requires std::ranges::range<V>;

        { v[n] } -> std::same_as<T&>;
        { *std::ranges::begin(v) } -> std::same_as<T&>;
    };

    template<typename V, typename T>
    concept vector2d = requires(V v, size_t n) {
        requires std::ranges::range<V>;

        { v[n] } -> vector1d<T>;
        { *std::ranges::begin(v) } -> vector1d<T>;
    };

    template<typename V, typename T>
    concept vector3d = requires(V v, size_t n) {
        requires std::ranges::range<V>;

        { v[n] } -> vector2d<T>;
        { *std::ranges::begin(v) } -> vector2d<T>;
    };

}// namespace ample::utils::concepts
