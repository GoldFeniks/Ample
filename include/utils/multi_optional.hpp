#pragma once
#include <tuple>
#include <optional>
#include <type_traits>

namespace ample::types {

    namespace _impl {

        template<typename, typename...>
        struct type_index;

        template<bool B, typename T, typename... Rest>
        struct type_index_helper {

            static constexpr size_t value = type_index<T, Rest...>::value + 1;

        };

        template<typename T, typename... Rest>
        struct type_index_helper<true, T, Rest...> {

            static constexpr size_t value = 0;

        };

        template<typename T, typename V, typename... Rest>
        struct type_index<T, V, Rest...> {

            static constexpr size_t value = type_index_helper<std::is_same_v<T, V>, T, Rest...>::value;

        };

        template<typename T>
        struct type_index<T> {

            static_assert(!std::is_same_v<T, T>, "Specified type was not found");

        };

        template<typename T, typename... Rest>
        static constexpr size_t type_index_v = type_index<T, Rest...>::value;

    }// namespace _impl

    template<typename... T>
    class multi_optional {

    public:

        template<typename... V>
        explicit constexpr multi_optional(V&&... values) {
            ((std::get<_impl::type_index_v<V, T...>>(_values) = std::forward<V>(values)), ...);
        }

        template<typename V>
        constexpr multi_optional& operator=(V&& value) {
            std::get<_impl::type_index_v<V, T...>>(_values) = std::forward<V>(value);
            return *this;
        }

        template<typename V>
        [[nodiscard]] constexpr bool has_value() const {
            return std::get<_impl::type_index_v<V, T...>>(_values).has_value();
        }

        template<typename V>
        constexpr const V& value() const {
            return std::get<_impl::type_index_v<V, T...>>(_values).value();
        }

    private:

        std::tuple<std::optional<T>...> _values;

    };

}// namespace ample::types

namespace std {

    template<typename V, typename... T>
    constexpr const V& get(const ample::types::multi_optional<T...>& optional) {
        return optional.template value<V>();
    }

}// namespace std