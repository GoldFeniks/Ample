#pragma once
#include <tuple>
#include <string>
#include <cstddef>
#include <fstream>
#include "modes.hpp"
#include "series.hpp"
#include "hydrology.hpp"
#include "bathymetry.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

namespace acstc {

    namespace __impl {

        template<typename T, typename V>
        struct modes_creator {

            static auto create(const json& data) {
                const auto type = data[0].template get<std::string>();
                if (type == "values") {
                    const auto x = data["/1/x"_json_pointer].template get<types::vector1d_t<T>>();
                    const auto y = data["/1/y"_json_pointer].template get<types::vector1d_t<T>>();
                    const auto& k_values = data["/1/k"_json_pointer];
                    types::vector3d_t<V> k_data(k_values.size(), types::vector2d_t<V>(x.size(), types::vector1d_t<V>(y.size())));
                    for (size_t i = 0; i < k_data.size(); ++i)
                        for (size_t j = 0; j < x.size(); ++j)
                            for (size_t k = 0; k < y.size(); ++k)
                                k_data[i][j][k] = V(k_values[i][j][2 * k], k_values[i][j][2 * k + 1]);
                    return std::make_tuple(
                            utils::linear_interpolated_data_2d<T, V>(x, y, std::move(k_data)),
                            utils::linear_interpolated_data_2d<T, T>(x, y,
                                    data["/1/phi"_json_pointer].template get<types::vector3d_t<T>>()));
                }
                if (type == "text_file")
                    return ::acstc::modes<T, V>::from_text(std::ifstream(data[1].template get<std::string>()));
                if (type == "binary_file")
                    return ::acstc::modes<T, V>::from_binary(std::ifstream(data[1].template get<std::string>(), std::ios::binary));
                throw std::logic_error("Unknown modes type: " + type);
            }

            static auto create_const(const json& data) {
                const auto type = data[0].template get<std::string>();
                if (type == "values") {
                    const auto y = data["/1/y"_json_pointer].template get<types::vector1d_t<T>>();
                    const auto& k_values = data["/1/k"_json_pointer];
                    types::vector2d_t<V> k_data(k_values.size(), types::vector1d_t<V>(y.size()));
                    for (size_t i = 0; i < k_data.size(); ++i)
                        for (size_t j = 0; j < y.size(); ++j)
                            k_data[i][j] = V(k_values[i][2 * j], k_values[i][2 * j + 1]);
                    return std::make_tuple(
                            utils::linear_interpolated_data_1d<T, V>(y, std::move(k_data)),
                            utils::linear_interpolated_data_1d<T, T>(y,
                                    data["/1/phi"_json_pointer].template get<types::vector2d_t<T>>()));
                }
                if (type == "text_file")
                    return ::acstc::modes<T, V>::const_from_text(std::ifstream(data[1].template get<std::string>()));
                if (type == "binary_file")
                    return ::acstc::modes<T, V>::const_from_binary(std::ifstream(data[1].template get<std::string>(), std::ios::binary));
                throw std::logic_error("Unknown modes type: " + type);
            }

        };

        template<typename T>
        struct modes_creator<T, T> {

            static auto create(const json& data) {
                const auto type = data[0].template get<std::string>();
                if (type == "values") {
                    const auto x = data["/1/x"_json_pointer].template get<types::vector1d_t<T>>();
                    const auto y = data["/1/y"_json_pointer].template get<types::vector1d_t<T>>();
                    return std::make_tuple(
                            utils::linear_interpolated_data_2d<T, T>(x, y,
                                    data["/1/k"_json_pointer].template get<types::vector3d_t<T>>()),
                            utils::linear_interpolated_data_2d<T, T>(x, y,
                                    data["/1/phi"_json_pointer].template get<types::vector3d_t<T>>()));
                }
                if (type == "text_file")
                    return ::acstc::modes<T, T>::from_text(std::ifstream(data[1].template get<std::string>()));
                if (type == "binary_file")
                    return ::acstc::modes<T, T>::from_binary(std::ifstream(data[1].template get<std::string>(), std::ios::binary));
                throw std::logic_error("Unknown modes type: " + type);
            }

            static auto create_const(const json& data) {
                const auto type = data[0].template get<std::string>();
                if (type == "values") {
                    const auto y = data["/1/y"_json_pointer].template get<types::vector1d_t<T>>();
                    return std::make_tuple(
                            utils::linear_interpolated_data_1d<T, T>(y,
                                    data["/1/k"_json_pointer].template get<types::vector2d_t<T>>()),
                            utils::linear_interpolated_data_1d<T, T>(y,
                                    data["/1/phi"_json_pointer].template get<types::vector2d_t<T>>()));
                }
                if (type == "text_file")
                    return ::acstc::modes<T, T>::const_from_text(std::ifstream(data[1].template get<std::string>()));
                if (type == "binary_file")
                    return ::acstc::modes<T, T>::const_from_binary(std::ifstream(data[1].template get<std::string>(), std::ios::binary));
                throw std::logic_error("Unknown modes type: " + type);
            }

        };

    }// namespace __impl;

#define CONFIG_DATA_FIELD(field, type) const type field () const { return _data[#field].template get<type>(); }
#define CONFIG_FIELD(field) const auto& field () const { return _##field; }

    template<typename T = types::real_t>
    class config {

    public:

        config() :
            _data(_default_data()),
            _bathymetry(_create_bathymetry(_data["bathymetry"])),
            _hydrology(_create_hydrology(_data["hydrology"]))
        {
            _fill_coefficients(_data["coefficients"]);
        }

        explicit config(const std::string& filename) : _data(_default_data()) {
            std::ifstream in(filename);
            json data;
            in >> data;
            _data.merge_patch(data);
            _bathymetry = _create_bathymetry(_data["bathymetry"]);
            _hydrology = _create_hydrology(_data["hydrology"]);
            _fill_coefficients(_data["coefficients"]);
        }

        CONFIG_DATA_FIELD(mode_subset, double)
        CONFIG_DATA_FIELD(x0, T)
        CONFIG_DATA_FIELD(x1, T)
        CONFIG_DATA_FIELD(nx, size_t)
        CONFIG_DATA_FIELD(y0, T)
        CONFIG_DATA_FIELD(y1, T)
        CONFIG_DATA_FIELD(ny, size_t)
        CONFIG_DATA_FIELD(ppm, size_t)
        CONFIG_DATA_FIELD(ordRich, size_t)
        CONFIG_DATA_FIELD(f, T)
        CONFIG_DATA_FIELD(z_r, T)
        CONFIG_DATA_FIELD(y_s, T)
        CONFIG_DATA_FIELD(complex_modes, bool)
        CONFIG_DATA_FIELD(const_modes, bool)
        CONFIG_DATA_FIELD(past_n, size_t)
        CONFIG_DATA_FIELD(border_width, size_t)
        CONFIG_DATA_FIELD(k0, types::vector1d_t<T>)
        CONFIG_DATA_FIELD(phi_s, types::vector1d_t<T>)

        CONFIG_FIELD(data)
        CONFIG_FIELD(a)
        CONFIG_FIELD(b)
        CONFIG_FIELD(c)

        auto x_bounds() const {
            return std::make_tuple(x0(), x1());
        }

        auto y_bounds() const {
            return std::make_tuple(y0(), y1());
        }

        auto bounds() const {
            return std::tuple_cat(x_bounds(), y_bounds());
        }

        const auto& bathymetry() const {
            return _bathymetry[0];
        }

        const auto& hydrology() const {
            return _hydrology[0];
        }

        auto coefficients() const {
            return std::make_tuple(a(), b(), c());
        }

        template<typename V = T>
        auto create_modes() const {
            if (_data.find("modes") != _data.end())
                return __impl::modes_creator<T, V>::create(_data["modes"]);
            return ::acstc::modes<T, V>::create(*this, 10, 10); //placeholder
        }

        template<typename V = T>
        auto create_const_modes() const {
            if (_data.find("modes") != _data.end())
                return __impl::modes_creator<T, V>::create_const(_data["modes"]);
            return ::acstc::modes<T, V>::create(*this, 10); //placeholder
        }

    private:

        json _data;
        T _a, _b, _c;
        utils::linear_interpolated_data_2d<T> _bathymetry;
        utils::delaunay_interpolated_data_2d<T> _hydrology;

        static json _default_data() {
            return {
                { "mode_subset", -1 },
                { "ppm", size_t(2) },
                { "ordRich", size_t(3) },
                { "f", T(100) },
                { "z_r", T(30) },
                { "y_s", T(0) },
                { "complex_modes", true },
                { "const_modes", true },
                { "past_n", size_t(0) },
                { "border_width", size_t(100) },
                { "bathymetry",
                  { "values",
                    {
                      { "x", { T(0), T(25000) } },
                      { "y", { T(-1000), T(1000) } },
                      { "values",
                        {
                          { T(25), T(375) },
                          { T(25), T(375) }
                        }
                      }
                    }
                  }
                },
                { "hydrology",
                  { "values",
                    {
                      { "x", { T(0), T(1) } },
                      { "z", { T(0), T(1) } },
                      { "values",
                        {
                          { T(1500), T(1500) },
                          { T(1500), T(1500) }
                        }
                      }
                    }
                  }
                },
                { "modes",
                  { "values",
                    {
                      { "y", { T(0), T(1) } },
                      { "k",
                        {
                          { T(0.1038243261),  T(7.31537e-6),   T(0.1038243261),  T(7.31537e-6)   },
                          { T(0.1010523078),  T(0.0000280559), T(0.1010523078),  T(0.0000280559) },
                          { T(0.09622526287), T(0.0000739264), T(0.09622526287), T(0.0000739264) },
                        }
                      },
                      { "phi",
                        {
                          { T(0.0372559), T(0.0372559) },
                          { T(0.0689824), T(0.0689824) },
                          { T(0.0882298), T(0.0882298) },
                        }
                      }
                    }
                  }
                },
                { "k0", { T(0.1038243261), T(0.1010523078), T(0.09622526287) } },
                { "phi_s", { T(0.0372559), T(0.0689824), T(0.0882298) } },
                { "x0", T(0) },
                { "x1", T(25000) },
                { "nx", size_t(2501) },
                { "y0", -T(1000) },
                { "y1",  T(1000) },
                { "ny", size_t(2001) },
                { "coefficients", { "pade" } },
            };
        }

        static auto _create_hydrology(const json& data) {
            const auto type = data[0].template get<std::string>();
            if (type == "values")
                return ::acstc::hydrology<T>::from_table(
                        data["/1/x"_json_pointer].template get<types::vector1d_t<T>>(),
                        data["/1/z"_json_pointer].template get<types::vector1d_t<T>>(),
                        data["/1/values"_json_pointer].template get<types::vector2d_t<T>>());
            if (type == "test_file")
                return ::acstc::hydrology<T>::from_text(std::ifstream(data[1].template get<std::string>()));
            if (type == "binary_file")
                return ::acstc::hydrology<T>::from_binary(std::ifstream(data[1].template get<std::string>(), std::ios::binary));
            throw std::logic_error("Unknown hydrology type: " + type);
        }

        static auto _create_bathymetry(const json& data) {
            const auto type = data[0].template get<std::string>();
            if (type == "values")
                return ::acstc::bathymetry<T>::from_table(
                        data["/1/x"_json_pointer].template get<types::vector1d_t<T>>(),
                        data["/1/y"_json_pointer].template get<types::vector1d_t<T>>(),
                        data["/1/values"_json_pointer].template get<types::vector2d_t<T>>());
            if (type == "text_file")
                return ::acstc::bathymetry<T>::from_text(std::ifstream(data[1].template get<std::string>()));
            if (type == "binary_file")
                return ::acstc::bathymetry<T>::from_binary(std::ifstream(data[1].template get<std::string>(), std::ios::binary));
            throw std::logic_error("Unknown bathymetry type: " + type);
        }

        void _fill_coefficients(const json& data) {
            const auto type = data[0].template get<std::string>();
            if (type == "pade") {
                const auto [a, b, c] = series::pade_series_coefficients<1, T>()[0];
                _a = a;
                _b = b;
                _c = c;
                return;
            }
            if (type == "abc") {
                _a = data["/1/a"_json_pointer].template get<T>();
                _b = data["/1/b"_json_pointer].template get<T>();
                _c = data["/1/c"_json_pointer].template get<T>();
                return;
            }
            throw std::logic_error("Unknown coefficients type: " + type);
        }

    };

}// namespace acstc
