#pragma once
#include <tuple>
#include <cstddef>
#include <type_traits>
#include "types.hpp"
#include "../io/reader.hpp"
#include "nlohmann/json.hpp"

namespace ample::utils {

    using nlohmann::json;

    template<typename T>
    struct var_dim;

    struct no_values_dim;

    namespace _impl {

        template<typename T>
        struct dim_vector {

            using type = types::vector1d_t<T>;

        };

        template<typename T>
        struct dim_vector<var_dim<T>> {

            using type = types::vector1d_t<typename dim_vector<T>::type>;

        };

        template<>
        struct dim_vector<no_values_dim> {

            using type = size_t;

        };

        template<typename T>
        using dim_vector_t = typename dim_vector<T>::type;

        template<typename D>
        struct size_getter {

            static size_t get(const types::vector1d_t<D>& data) {
                return data.size();
            }

        };

        template<typename D>
        struct size_getter<var_dim<D>> {

            template<typename... Args>
            static size_t get(const dim_vector_t<var_dim<D>>& data, const size_t& index, const Args&... args) {
                return size_getter<D>::get(data[index], args...);
            }

            static size_t get(const dim_vector_t<var_dim<D>>& data) {
                return data.size();
            }

        };

        template<>
        struct size_getter<no_values_dim> {

            static size_t get(const size_t& size) {
                return size;
            }

        };

        template<typename D>
        struct value_getter {

            static const auto& get(const types::vector1d_t<D>& data) {
                return data;
            }

        };

        template<typename D>
        struct value_getter<var_dim<D>> {

            template<typename... Args>
            static const auto& get(const dim_vector_t<var_dim<D>>& data, const size_t& index, const Args&... args) {
                return value_getter<D>::get(data[index], args...);
            }

            static const auto& get(const dim_vector_t<var_dim<D>>& data) {
                return data;
            }

        };

        template<>
        struct value_getter<no_values_dim> {

            template<typename D>
            static auto get(const D&, ...) {
                static_assert(!std::is_same_v<D, D>, "Cannot get values for no values dimension");
            }

        };

        template<typename D>
        struct data_reader {

            static auto read(const json& data) {

                using vector_type = types::vector1d_t<D>;
                using mesh_type = decltype(utils::mesh_1d<D>(std::declval<D>(), std::declval<D>(), std::declval<size_t>()));

                if (data.is_number())
                    return utils::mesh_1d(D(0), D(1), data.template get<size_t>());

                if (data.contains("values") && data.contains("bounds"))
                    throw std::logic_error("Found both values and bounds keys");

                const auto n = data["n"].template get<size_t>();

                if (data.contains("values")) {
                    if (data["values"].is_string()) {
                        std::ifstream inp(data["values"].get<std::string>());
                        return make_vector_i(n, [&inp](auto&&) {
                                D v;
                                inp >> v;
                                return v;
                            }
                        );
                    }

                    const auto result = data["values"].get<vector_type>();
                    if (result.size() != n)
                        throw std::logic_error(
                            std::string("Wrong values specification size. Got ") +
                            std::to_string(result.size()) +
                            ". Required " +
                            std::to_string(n)
                        );
                    return result;
                }

                if constexpr (std::is_same_v<vector_type, mesh_type>) {
                    if (data.contains("bounds"))
                        return utils::mesh_1d(
                            data["/bounds/a"_json_pointer].get<D>(),
                            data["/bounds/b"_json_pointer].get<D>(),
                            n
                        );
                } else if (data.contains("bounds"))
                    throw std::logic_error("Bounds specification is unsupported for underlying type");

                return utils::mesh_1d(D(0), D(1), n);
            }

        };

        template<typename D>
        struct data_reader<var_dim<D>> {

            static auto read(const json& data) {
                if (data.is_array())
                    return utils::make_vector(data, data_reader<D>::read);
                return dim_vector_t<var_dim<D>>{ data_reader<D>::read(data) };
            }
        };

        template<>
        struct data_reader<no_values_dim> {

            static auto read(const json& data) {
                if (data.is_object()) {
                    if (data.contains("values") || data.contains("bounds"))
                        throw std::logic_error("Cannot specify values for no values dimension");
                    return data["n"].template get<size_t>();
                }

                if (data.is_number())
                    return data.template get<size_t>();

                throw std::logic_error("No value dimension specification must be either a single number or an object with \"n\" key");
            }

        };

    }// namespace _impl

    template<typename>
    struct is_variable_dim {

        static constexpr bool value = false;

    };

    template<typename T>
    struct is_variable_dim<var_dim<T>> {

        static constexpr bool value = true;

    };

    template<typename T>
    constexpr bool is_variable_dim_v = is_variable_dim<T>::value;

    template<typename... T>
    class dimensions {

    private:

        using value_t = std::tuple<T...>;
        using data_t = std::tuple<_impl::dim_vector_t<T>...>;

    public:

        static constexpr size_t N = sizeof...(T);

        explicit dimensions(const json& data) {
            if (data.size() != N)
                throw std::runtime_error(
                    std::string("Wrong number of dimensions. Got ") +
                    std::to_string(data.size()) +
                    ". Required " +
                    std::to_string(N)
                );

            _read_dimensions(data, std::make_index_sequence<N>{});
        }

        template<size_t N, typename... Args>
        const auto& get(const Args&... args) const {
            return _impl::value_getter<std::tuple_element_t<N, value_t>>::get(std::get<N>(_data), args...);
        }

        template<size_t N, typename... Args>
        auto size(const Args&... args) const {
            return _impl::size_getter<std::tuple_element_t<N, value_t>>::get(std::get<N>(_data), args...);
        }

        template<size_t M>
        static constexpr bool is_variable_dim = is_variable_dim_v<std::tuple_element_t<M, value_t>>;

    private:

        data_t _data;

        template<size_t... I>
        void _read_dimensions(const json& data, const std::integer_sequence<size_t, I...>&) {
            ((std::get<I>(_data) = _impl::data_reader<std::tuple_element_t<I, value_t>>::read(data[I])),...);
        }

    };

}// namespace ample::utils
