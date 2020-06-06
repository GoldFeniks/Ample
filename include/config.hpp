#pragma once
#include <cmath>
#include <tuple>
#include <string>
#include <cstddef>
#include <fstream>
#include <filesystem>
#include <unordered_map>
#include "modes.hpp"
#include "series.hpp"
#include "hydrology.hpp"
#include "bathymetry.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "nlohmann/json.hpp"
#include "initial_conditions.hpp"

namespace acstc {

    using nlohmann::json;

    namespace __impl {

        template<typename T, typename V>
        struct modes_creator {

            static auto create(const json& data, const size_t count, const std::filesystem::path& path) {
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
                    auto phi_data = data["/1/phi"_json_pointer].template get<types::vector3d_t<T>>();
                    ::acstc::modes<T, V>::smooth((y.back() - y.front()) / (y.size() - 1), count, k_data, phi_data);
                    return std::make_tuple(
                            utils::linear_interpolated_data_2d<T, V>(x, y, std::move(k_data)),
                            utils::linear_interpolated_data_2d<T, T>(x, y, std::move(phi_data)));
                }
                if (type == "text_file")
                    return ::acstc::modes<T, V>::from_text(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>())), count);
                if (type == "binary_file")
                    return ::acstc::modes<T, V>::from_binary(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>()), std::ios::binary), count);
                throw std::logic_error("Unknown modes type: " + type);
            }

            static auto create_const(const json& data, const std::filesystem::path& path) {
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
                    return ::acstc::modes<T, V>::const_from_text(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>())));
                if (type == "binary_file")
                    return ::acstc::modes<T, V>::const_from_binary(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>()), std::ios::binary));
                throw std::logic_error("Unknown modes type: " + type);
            }

        };

        template<typename T>
        struct modes_creator<T, T> {

            static auto create(const json& data, const size_t count, const std::filesystem::path& path) {
                const auto type = data[0].template get<std::string>();
                if (type == "values") {
                    const auto x = data["/1/x"_json_pointer].template get<types::vector1d_t<T>>();
                    const auto y = data["/1/y"_json_pointer].template get<types::vector1d_t<T>>();
                    auto k_data = data["/1/k"_json_pointer].template get<types::vector3d_t<T>>();
                    auto phi_data = data["/1/phi"_json_pointer].template get<types::vector3d_t<T>>();
                    ::acstc::modes<T, T>::smooth((y.back() - y.front()) / (y.size() - 1), count, k_data, phi_data);
                    return std::make_tuple(
                            utils::linear_interpolated_data_2d<T, T>(x, y, std::move(k_data)),
                            utils::linear_interpolated_data_2d<T, T>(x, y, std::move(phi_data)));
                }
                if (type == "text_file")
                    return ::acstc::modes<T, T>::from_text(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>())), count);
                if (type == "binary_file")
                    return ::acstc::modes<T, T>::from_binary(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>()), std::ios::binary), count);
                throw std::logic_error("Unknown modes type: " + type);
            }

            static auto create_const(const json& data, const std::filesystem::path& path) {
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
                    return ::acstc::modes<T, T>::const_from_text(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>())));
                if (type == "binary_file")
                    return ::acstc::modes<T, T>::const_from_binary(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>()), std::ios::binary));
                throw std::logic_error("Unknown modes type: " + type);
            }

        };

    }// namespace __impl;

#define CONFIG_DATA_FIELD(field, type)                                                   \
    const type& field () const {                                                         \
        if (!_cache.count(#field))                                                       \
            _cache[#field] = new data_field<type>(_data[#field].template get<type>());   \
        return _cache[#field]->template cast<type>().value;                              \
    }                                                                                    \
    void field(const type& value) {                                                      \
        if (const auto it = _cache.find(#field); it != _cache.end())                     \
            _cache.erase(it);                                                            \
        _data[#field] = value;                                                           \
    }                                                                                    \
    void field(type&& value) {                                                           \
        if (const auto it = _cache.find(#field); it != _cache.end())                     \
            _cache.erase(it);                                                            \
        _data[#field] = std::move(value);                                                \
    }


#define CONFIG_FIELD(field, type)                                                        \
    const type& field () const { return _##field; }                                      \
    void field (const type& value) { _##field = value; }                                 \
    void field (type&& value) { _##field = std::move(value); }


    template<typename T = types::real_t>
    class config {

    public:

        config() :
            _data(_default_data()),
            _bathymetry(_create_bathymetry(_data["bathymetry"], "")),
            _hydrology(_create_hydrology(_data["hydrology"], ""))
        {
            _fill_coefficients(_data["coefficients"]);
        }

        explicit config(const std::string& filename) : _data(_default_data()) {
            update_from_file(filename);
        }

        ~config() {
            for (const auto& [key, value] : _cache)
                delete value;
        }

        void update_from_file(const std::string& filename) {
            _path = filename;
            std::ifstream in(filename);
            json data;
            in >> data;
            _data.merge_patch(data);
            _bathymetry = _create_bathymetry(_data["bathymetry"], _path);
            _hydrology = _create_hydrology(_data["hydrology"], _path);
            _fill_coefficients(_data["coefficients"]);

        }

        CONFIG_DATA_FIELD(mode_subset, double)
        CONFIG_DATA_FIELD(max_mode, size_t)
        CONFIG_DATA_FIELD(n_modes, size_t)
        CONFIG_DATA_FIELD(x0, T)
        CONFIG_DATA_FIELD(x1, T)
        CONFIG_DATA_FIELD(nx, size_t)
        CONFIG_DATA_FIELD(y0, T)
        CONFIG_DATA_FIELD(y1, T)
        CONFIG_DATA_FIELD(ny, size_t)
        CONFIG_DATA_FIELD(ppm, size_t)
        CONFIG_DATA_FIELD(ordRich, size_t)
        CONFIG_DATA_FIELD(f, T)
        CONFIG_DATA_FIELD(z_s, T)
        CONFIG_DATA_FIELD(z_r, T)
        CONFIG_DATA_FIELD(y_s, T)
        CONFIG_DATA_FIELD(n_layers, size_t)
        CONFIG_DATA_FIELD(bottom_rhos, types::vector1d_t<T>)
        CONFIG_DATA_FIELD(betas, types::vector1d_t<T>)
        CONFIG_DATA_FIELD(bottom_layers, types::vector1d_t<T>)
        CONFIG_DATA_FIELD(bottom_c1s, types::vector1d_t<T>)
        CONFIG_DATA_FIELD(bottom_c2s, types::vector1d_t<T>)
        CONFIG_DATA_FIELD(complex_modes, bool)
        CONFIG_DATA_FIELD(const_modes, bool)
        CONFIG_DATA_FIELD(additive_depth, bool)
        CONFIG_DATA_FIELD(past_n, size_t)
        CONFIG_DATA_FIELD(border_width, size_t)
        CONFIG_DATA_FIELD(k0, types::vector1d_t<T>)
        CONFIG_DATA_FIELD(phi_s, types::vector1d_t<T>)
        CONFIG_DATA_FIELD(a0, T)
        CONFIG_DATA_FIELD(a1, T)
        CONFIG_DATA_FIELD(na, size_t)
        CONFIG_DATA_FIELD(l0, T)
        CONFIG_DATA_FIELD(l1, T)
        CONFIG_DATA_FIELD(nl, size_t)
        CONFIG_DATA_FIELD(init, std::string)

        CONFIG_FIELD(data, json)
        CONFIG_FIELD(a, T)
        CONFIG_FIELD(b, T)
        CONFIG_FIELD(c, T)

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
        auto create_modes(const size_t& c = 0, const bool show_progress = false) const {
            if (_data.count("modes"))
                return __impl::modes_creator<T, V>::create(_data["modes"], border_width(), _path);
            if (_data.count("mnx") && _data.count("mny"))
                return ::acstc::modes<T, V>::create(*this, _data["mnx"].template get<size_t>(), _data["mny"].template get<size_t>(), c, show_progress);
            return ::acstc::modes<T, V>::create(*this, c, show_progress);
        }

        template<typename V = T>
        auto create_const_modes(const size_t& c = 0, const bool show_progress = false) const {
            if (_data.count("modes"))
                return __impl::modes_creator<T, V>::create_const(_data["modes"], _path);
            if (_data.count("mny"))
                return ::acstc::modes<T, V>::create(*this, _data["mny"].template get<size_t>(), c, show_progress);
            return ::acstc::modes<T, V>::create(*this, bathymetry().y().size(), c, show_progress);
        }

        auto create_source_modes(const size_t& c = 0) const {
            if (_data.count("k0") && _data.count("phi_s"))
                return std::make_tuple(k0(), phi_s());
            const auto x = T(0);
            auto n_m = ::acstc::modes<T>::calc_modes(*this, x, bathymetry().point(x, y_s()), z_s(), c);
            types::vector1d_t<T> k0(n_m.khs.size()), phi_s(n_m.khs.size());
            for (size_t i = 0; i < k0.size(); ++i) {
                k0[i] = n_m.khs[i];
                phi_s[i] = n_m.mfunctions_zr[i][0];
            }
            return std::make_tuple(k0, phi_s);
        }

        auto get_tapering_parameters() const {
            const auto type = _data["/tapering/0"_json_pointer].get<std::string>();
            const auto& data = _data["/tapering/1"_json_pointer];
            if (data.count("left") && data.count("right"))
                return std::make_tuple(type, data["left"].get<T>(), data["right"].get<T>());
            const auto value = data["value"].get<T>();
            return std::make_tuple(type, value, value);
        }

    private:

        template<typename>
        class data_field;

        struct data_field_base {

            template<typename C>
            const data_field<C>& cast() const {
                return *static_cast<const data_field<C>*>(this);
            }

        };

        template<typename C>
        class data_field : public data_field_base {

        public:

            const C value;

            data_field(const C& value) : value(value) {}
            data_field(C&& value) : value(std::move(value)) {}

        };

        json _data;
        T _a, _b, _c;

        std::filesystem::path _path;
        utils::linear_interpolated_data_2d<T> _bathymetry;
        utils::delaunay_interpolated_data_2d<T> _hydrology;
        mutable std::unordered_map<std::string, data_field_base*> _cache;

        static json _default_data() {
            return {
                { "mode_subset", -1 },
                { "max_mode", size_t(-1) },
                { "n_modes", size_t(0) },
                { "ppm", size_t(2) },
                { "ordRich", size_t(3) },
                { "f", T(25) },
                { "z_r", T(30) },
                { "z_s", T(100) },
                { "y_s", T(0) },
                { "n_layers", size_t(1) },
                { "bottom_layers", { T(500) } },
                { "bottom_c1s", { T(1700) } },
                { "bottom_c2s", { T(1700) } },
                { "bottom_rhos", { T(1.5) } },
                { "betas", { T(0), T(0.5) } },
                { "complex_modes", true },
                { "const_modes", true },
                { "additive_depth", false },
                { "past_n", size_t(0) },
                { "border_width", size_t(10) },
                { "bathymetry",
                  { "values",
                    {
                      { "x", { T(0), T(1) } },
                      { "y", { T(0), T(1) } },
                      { "values",
                        {
                          { T(200), T(200) },
                          { T(200), T(200) }
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
                { "x0", T(0) },
                { "x1", T(15000) },
                { "nx", size_t(15001) },
                { "y0", -T(4000) },
                { "y1",  T(4000) },
                { "ny", size_t(8001) },
                { "coefficients", { "pade" } },
                { "a0", -T(M_PI) / T(4) },
                { "a1",  T(M_PI) / T(4) },
                { "na", size_t(90) },
                { "l0", T(0) },
                { "l1", T(4000) },
                { "nl", size_t(4001) },
                { "init", "green" },
                { "tapering", 
                  { "angled",
                    {
                      { "value", T(0.1) }
                    }
                  }
                }
            };
        }

        static auto _create_hydrology(const json& data, const std::filesystem::path& path) {
            const auto type = data[0].template get<std::string>();
            if (type == "values")
                return ::acstc::hydrology<T>::from_table(
                        data["/1/x"_json_pointer].template get<types::vector1d_t<T>>(),
                        data["/1/z"_json_pointer].template get<types::vector1d_t<T>>(),
                        data["/1/values"_json_pointer].template get<types::vector2d_t<T>>());
            if (type == "text_file")
                return ::acstc::hydrology<T>::from_text(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>())));
            if (type == "binary_file")
                return ::acstc::hydrology<T>::from_binary(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>()), std::ios::binary));
            throw std::logic_error("Unknown hydrology type: " + type);
        }

        static auto _create_bathymetry(const json& data, const std::filesystem::path& path) {
            const auto type = data[0].template get<std::string>();
            if (type == "values")
                return ::acstc::bathymetry<T>::from_table(
                        data["/1/x"_json_pointer].template get<types::vector1d_t<T>>(),
                        data["/1/y"_json_pointer].template get<types::vector1d_t<T>>(),
                        data["/1/values"_json_pointer].template get<types::vector2d_t<T>>());
            if (type == "text_file")
                return ::acstc::bathymetry<T>::from_text(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>())));
            if (type == "binary_file")
                return ::acstc::bathymetry<T>::from_binary(std::ifstream(utils::make_file_path(path, data[1].template get<std::string>()), std::ios::binary));
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
