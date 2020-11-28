#pragma once
#include <cmath>
#include <tuple>
#include <string>
#include <cstddef>
#include <fstream>
#include <utility>
#include <optional>
#include <filesystem>
#include <type_traits>
#include <unordered_map>
#include "modes.hpp"
#include "series.hpp"
#include "hydrology.hpp"
#include "bathymetry.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "nlohmann/json.hpp"
#include "utils/convertors.hpp"
#include "utils/dimensions.hpp"
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

        template<typename T>
        auto create_receiver_depth(const types::vector1d_t<std::array<T, 3>>& receivers) {
            types::vector1d_t<T> values;
            types::vector1d_t<std::tuple<T, T>> points;
            for (const auto& [x, y, z] : receivers) {
                values.emplace_back(z);
                points.emplace_back(x, y);
            }
            return utils::nearest_neighbour_interpolated_data_2d<T>(std::move(points), std::move(values));
        }

        template<typename T>
        auto create_range(const json& data) {
            return acstc::utils::mesh_1d(
                data["a"].template get<T>(), 
                data["b"].template get<T>(),
                data["n"].template get<size_t>());
        }

        template<typename T, typename... D>
        struct input_data {

            utils::dimensions<D...> dimensions;
            decltype(acstc::vector_reader<T>::template read<0, D...>(std::declval<std::istream>(), dimensions)) data;

            explicit input_data(const json& data, const std::filesystem::path& path) :
                dimensions(data["dimensions"]),
                data(get_data(data["values"], dimensions, path, data.contains("binary") && data["binary"].template get<bool>())) {}

            template<size_t M = 0>
            static auto get_data(const json& data, const utils::dimensions<D...>& dims, 
                const std::filesystem::path& path, const bool& binary) {

                if (data.is_string())
                    return read_data<M>(dims, path, data.template get<std::string>(), binary);

                check_size<M>(data.size(), dims.template size<M>());

                if constexpr (M + 1 < sizeof...(D))
                    if constexpr (utils::dimensions<D...>::template is_variable_dim<M>)
                        return utils::make_vector_i(data,
                            [&dims, &path, &binary](const auto& data, const size_t& i) {
                                check_size<M>(data.size(), dims.template size<M>(i));
                                return make_vector<M>(data, dims, path, binary);
                            }
                        );
                    else 
                        return make_vector<M>(data, dims, path, binary);
                else
                    if constexpr (utils::dimensions<D...>::template is_variable_dim<M>)
                        return utils::make_vector_i(data,
                            [&dims, &path, &binary](const auto& data, const size_t& i) {
                                if (data.is_string()) {
                                    types::vector1d_t<T> result(dims.template size<M>(i));
                                    if (binary) {
                                        std::ifstream inp(utils::make_file_path(path, data.template get<std::string>()), std::ios::binary);
                                        inp.read(reinterpret_cast<char*>(result.data()), sizeof(T) * result.size());
                                    } else {
                                        std::ifstream inp(utils::make_file_path(path, data.template get<std::string>()));
                                        for (auto& it : result)
                                            inp >> it;
                                    }
                                    return result;
                                }

                                check_size<M>(data.size(), dims.template size<M>(i));
                                return data.template get<types::vector1d_t<T>>();
                            }
                        );
                    else 
                        return data.template get<types::vector1d_t<T>>();
            }

            template<size_t M>
            static auto check_size(const size_t& n, const size_t& m) {
                if (n != m)
                    throw std::runtime_error(
                        std::string("Wrong dimension size at level ") +
                        std::to_string(M) + 
                        ". Expected " +
                        std::to_string(m) + 
                        ". Got " +
                        std::to_string(n)
                    );
            }

            template<size_t M>
            static auto read_data(const utils::dimensions<D...>& dims, const std::filesystem::path& path, const std::string& filename, const bool& binary) {
                const auto file_path = utils::make_file_path(path, filename);
                if (binary)
                    return acstc::vector_reader<T>::template binary_read<M, D...>(std::ifstream(file_path, std::ios::binary), dims);
                return acstc::vector_reader<T>::template read<M, D...>(std::ifstream(file_path), dims);
            }

            template<size_t M>
            static auto make_vector(const json& data, const utils::dimensions<D...>& dims, 
                const std::filesystem::path& path, const bool& binary) {
                return utils::make_vector(data, 
                    [&dims, &path, &binary](const auto& data) {
                        return get_data<M + 1>(data, dims, path, binary);
                    }
                );
            }

        };

    }// namespace __impl;

#define CONFIG_DATA_FIELD(field, type)                                      \
    private:                                                                \
        mutable std::optional<type> _data_##field;                          \
    public:                                                                 \
        const type& field () const {                                        \
            if (!_data_##field.has_value())                                 \
                _data_##field = type(_data[#field].template get<type>());   \
            return _data_##field.value();                                   \
        }                                                                   \
        void field(const type& value) {                                     \
            _data_##field = value;                                          \
        }                                                                   \
        void field(type&& value) {                                          \
            _data_##field = std::move(value);                               \
        }

#define CONFIG_FIELD(field, type)                                           \
    private:                                                                \
        type _##field;                                                      \
    public:                                                                 \
        const type& field () const { return _##field; }                     \
        void field (const type& value) { _##field = value; }                \
        void field (type&& value) { _##field = std::move(value); }

#define CONFIG_INPUT_DATA(field, type, op)                           \
    private:                                                                \
        std::optional<type> _##field;                                       \
    public:                                                                 \
        const auto& field () const {                                        \
            if (!_##field.has_value())                                      \
                throw std::logic_error("Field " #field " has no value");    \
            return _##field.value() op;                                     \
        }                                                                   \
        bool has_##field () const {                                         \
            return _##field.has_value();                                    \
        }                                                                   \
        template<typename... Args>                                          \
        std::enable_if_t<0 < sizeof...(Args)> field(Args&&... args) {       \
            _##field = type(std::forward<Args>(args)...);                   \
        }


    template<typename T = types::real_t>
    class config {

    public:

        config() :
            _data(_default_data())
        {
            _fill_coefficients(_data["coefficients"]);
        }

        explicit config(const std::string& filename) : _data(_default_data()) {
            update_from_file(filename);
        }

        void update_from_file(const std::string& filename) {
            _path = filename;
            std::ifstream in(filename);
            json data;
            in >> data;
            _data.merge_patch(data);

            for (const auto& it : _data["input_data"]) {
                const auto type = it["type"].template get<std::string>();

                try {
                    if (type == "bathymetry") {
                        _bathymetry = _create_bathymetry(it, _path);
                        continue;
                    }

                    if (type == "hydrology") {
                        _hydrology = _create_hydrology(it, _path);
                        continue;
                    }

                    if (type == "receivers") {
                        _receiver_depth = _create_receiver_depth(it, _path);
                        continue;
                    }

                    if (type == "source_function") {
                        std::tie(_times, _source_function) = _create_source_function<T>(it, _path);
                        continue;
                    }

                    if (type == "source_spectrum") {
                        if (_frequencies.has_value())
                            throw std::logic_error("Multiple values given for frequencies");

                        std::tie(_frequencies, _source_spectrum) = _create_source_function<std::complex<T>>(it, _path);
                        continue;
                    }

                    if (type == "frequencies") {
                        if (_frequencies.has_value())
                            throw std::logic_error("Multiple values given for frequencies");

                        _frequencies = _create_frequences(it, _path);
                        continue;
                    }
                } catch (const std::runtime_error& error) {
                    throw std::runtime_error(utils::join("Error while parsing ", type, ": ", error.what()));
                }

                throw std::runtime_error(std::string("Unknown input data type \"") + type + "\"");
            }
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
        CONFIG_DATA_FIELD(ord_rich, size_t)
        CONFIG_DATA_FIELD(z_s, T)
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
        CONFIG_DATA_FIELD(tolerance, T)
        CONFIG_DATA_FIELD(reference_index, size_t)
        CONFIG_DATA_FIELD(sel_range, types::tuple2_t<T>)
        CONFIG_DATA_FIELD(sel_strict, bool)

        CONFIG_FIELD(data, json)
        CONFIG_FIELD(a, T)
        CONFIG_FIELD(b, T)
        CONFIG_FIELD(c, T)

        CONFIG_INPUT_DATA(bathymetry, utils::linear_interpolated_data_2d<T>, [0])
        CONFIG_INPUT_DATA(hydrology, utils::delaunay_interpolated_data_2d<T>, [0])
        CONFIG_INPUT_DATA(receiver_depth, utils::nearest_neighbour_interpolated_data_2d<T>, [0])
        CONFIG_INPUT_DATA(source_function, types::vector1d_t<T>, ;)
        CONFIG_INPUT_DATA(source_spectrum, types::vector1d_t<std::complex<T>>, ;)
        CONFIG_INPUT_DATA(frequencies, types::vector1d_t<T>, ;)
        CONFIG_INPUT_DATA(times, types::vector1d_t<T>, ;)

        auto x_bounds() const {
            return std::make_tuple(x0(), x1());
        }

        auto y_bounds() const {
            return std::make_tuple(y0(), y1());
        }

        auto bounds() const {
            return std::tuple_cat(x_bounds(), y_bounds());
        }

        auto coefficients() const {
            return std::make_tuple(a(), b(), c());
        }

        auto f() const {
            return frequencies()[_index];
        }

        auto t() const {
            return times()[_index];
        }

        auto dt() const {
            return times().size() >= 2 ? times()[1] - times()[0] : T(0);
        }

        void index(const size_t& index) {
            _index = index;
        }

        const auto& index() const {
            return _index;
        }

        auto z_r(const T& x = T(0), const T& y = T(0)) const {
            return receiver_depth()[0].point(x, y);
        }

        auto mnx() const {
            if (_data.count("mnx") && _data.count("mny"))
                return _data["mnx"].template get<size_t>();
            return bathymetry().x().size();
        }

        auto mny() const {
            if (_data.count("mny") && (const_modes() || _data.count("mnx")))
                return _data["mny"].template get<size_t>();
            return bathymetry().y().size();
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

        void save(const std::filesystem::path& output) const {
            json out = _data;

            for (auto& it : out["input_data"])
                _copy_files(it["values"], _get_dim_count(it["dimensions"]), _path, 
                    output / it["type"].template get<std::string>(), it.contains("binary") && it["binary"].template get<bool>());

            if (!_data.count("mnx") || !_data.count("mny")) {
                out["mnx"] = bathymetry().x().size();
                out["mny"] = bathymetry().y().size();
            }

            const auto& coeffs = _data["coefficients"];
            if (!(coeffs[0].template get<std::string>() == "abc"))
                out["coefficients"] = { "abc", { a(), b(), c() } };

            std::ofstream file(output / "config.json");
            file << std::setw(4) << out << std::endl;
        }

    private:

        mutable size_t _index = 0;

        std::filesystem::path _path;

        static json _default_data() {
            return {
                { "mode_subset", -1 },
                { "max_mode", -1 },
                { "n_modes", size_t(0) },
                { "ppm", size_t(2) },
                { "ord_rich", size_t(3) },
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
                { "init", "greene" },
                { "tapering", 
                  { "angled",
                    {
                      { "value", T(0.1) }
                    }
                  }
                },
                { "tolerance", T(0.02) },
                { "reference_index", size_t(0) },
                { "sel_range", { T(-1), T(-1) } },
                { "sel_strict", false }
            };
        }

        static size_t _get_depth(const json& data) {
            return 1 + (data.is_array() ? _get_depth(data[0]) : 0);
        }

        static size_t _get_dim_count(const json& data) {
            size_t result = 0;
            for (const auto& it : data)
                result += _get_depth(it);
            return result;
        }

        static void _copy_files(json& data, const size_t& depth, 
            const std::filesystem::path& path, const std::filesystem::path& output, const bool& binary) {
            size_t count = 0;
            _copy_files(data, count, depth, path, output, binary);
        }

        static void _copy_files(json& data, size_t& count, const size_t& depth,
            const std::filesystem::path& path, const std::filesystem::path& output, const bool& binary) {
            if (depth == 0)
                return;

            if (data.is_string()) {
                std::filesystem::create_directories(output);

                auto filename = output / std::to_string(count++);
                filename += binary ? ".bin" : ".txt";
                std::filesystem::copy(utils::make_file_path(path, data.template get<std::string>()), filename,
                    std::filesystem::copy_options::overwrite_existing);

                data = filename;
                return;
            }

            if (data.is_array()) {
                for (auto& it : data)
                    _copy_files(it, count, depth - 1, path, output, binary);

                return;
            }

            throw std::logic_error("Data must be either string or array");
        }

        static auto _create_hydrology(const json& data, const std::filesystem::path& path) {
            const auto [dimensions, inp_data] = __impl::input_data<T, T, T>(data, path);

            return ::acstc::hydrology<T>::from_table(
                dimensions.template get<1>(),
                dimensions.template get<0>(),
                inp_data
            );
        }

        static auto _create_bathymetry(const json& data, const std::filesystem::path& path) {
            const auto [dimensions, inp_data] = __impl::input_data<T, T, T>(data, path);

            return ::acstc::bathymetry<T>::from_table(
                dimensions.template get<0>(),
                dimensions.template get<1>(),
                inp_data
            );
        }

        static auto _create_receiver_depth(const json& data, const std::filesystem::path& path) {
            return __impl::create_receiver_depth(__impl::input_data<std::array<T, 3>, utils::no_values_dim>(data, path).data);
        }

        template<typename V = T>
        static auto _create_source_function(const json& data, const std::filesystem::path& path) {
            const auto [dimensions, inp_data] = __impl::input_data<V, T>(data, path);
            return std::make_tuple(dimensions.template get<0>(), inp_data);
        }

        static auto _create_frequences(const json& data, const std::filesystem::path& path) {
            return __impl::input_data<T, utils::no_values_dim>(data, path).data;
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
