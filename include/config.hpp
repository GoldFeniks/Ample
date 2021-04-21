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
#include "io/reader.hpp"
#include "feniks/zip.hpp"
#include "utils/join.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "nlohmann/json.hpp"
#include "utils/convertors.hpp"
#include "utils/dimensions.hpp"
#include "initial_conditions.hpp"
#include "utils/multi_optional.hpp"

namespace ample {

    using nlohmann::json;

    namespace _impl {

        HAS_CONCEPT(can_make_mesh, utils::mesh_1d(std::declval<C>(), std::declval<C>(), size_t(0)), typename = void)

        template<typename T, typename... D>
        struct input_data {

            utils::dimensions<D...> dimensions;
            decltype(ample::vector_reader<T>::template read<0, D...>(std::declval<std::istream>(), dimensions)) data;

            explicit input_data(const json& data, const std::filesystem::path& path) :
                dimensions(data["dimensions"]),
                data(get_data(data["values"], dimensions, path, data.contains("binary") && data["binary"].template get<bool>())) {}

            template<size_t M = 0>
            static auto get_data(const json& data, const utils::dimensions<D...>& dims, 
                const std::filesystem::path& path, const bool& binary) {

                if (data.is_string())
                    return read_data<M>(dims, path, data.template get<std::string>(), binary);

                if constexpr (M + 1 < sizeof...(D)) {
                    check_array_size<M>(data, dims.template size<M>());

                    if constexpr (utils::dimensions<D...>::template is_variable_dim<M>)
                        return utils::make_vector_i(data,
                            [&dims, &path, &binary](const auto &data, const size_t &i) {
                                check_size<M>(data.size(), dims.template size<M>(i));
                                return make_vector<M>(data, dims, path, binary);
                            }
                        );
                    else
                        return make_vector<M>(data, dims, path, binary);
                } else {
                    if constexpr (utils::dimensions<D...>::template is_variable_dim<M>) {
                        check_array_size<M>(data, dims.template size<M>());

                        return utils::make_vector_i(data,
                            [&dims, &path, &binary](const auto& data, const size_t& i) {
                                if (data.is_string()) {
                                    if (binary) {
                                        types::vector1d_t<T> result(dims.template size<M>(i));
                                        std::ifstream inp(utils::make_file_path(path, data.template get<std::string>()), std::ios::binary);

                                        inp.read(reinterpret_cast<char*>(result.data()), sizeof(T) * result.size());
                                        utils::dynamic_assert(inp.gcount() == sizeof(T) * dims.template size<M>(i), "Insufficient data in file");
                                        return result;
                                    }

                                    types::vector1d_t<T> result;
                                    std::ifstream inp(utils::make_file_path(path, data.template get<std::string>()));

                                    _impl::read_line(inp, result);
                                    utils::dynamic_assert(result.size() == dims.template size<M>(i),
                                                          "Incorrect number of elements on line 1. Expected ",
                                                          dims.template size<M>(i), ", but got ", result.size());
                                    return result;
                                }

                                if constexpr (can_make_mesh_v<T>)
                                    if (data.is_object()) {
                                        utils::dynamic_assert(data.contains("a") && data.contains("b") && data.size() == 2 &&
                                                              data["a"].is_number() && data["b"].is_number(),
                                                              "Expected ranged data specification as { \"a\": <number>, \"b\": <number> }, but got ", data, " at level ", M);

                                        return utils::mesh_1d(data["a"].template get<T>(), data["b"].template get<T>(), dims.template size<M>(i));
                                    }

                                check_array_size<M>(data, dims.template size<M>(i));
                                return data.template get<types::vector1d_t<T>>();
                            }
                        );
                    } else {
                        if constexpr (can_make_mesh_v<T>)
                            if (data.is_object()) {
                                utils::dynamic_assert(data.contains("a") && data.contains("b") && data.size() == 2,
                                                      data["a"].is_number() && data["b"].is_number(),
                                                      "Expected ranged data specification as { \"a\": <number>, \"b\": <number> }, but got ", data, " at level ", M);

                                return utils::mesh_1d(data["a"].template get<T>(), data["b"].template get<T>(), dims.template size<M>());
                            }

                        check_array_size<M>(data, dims.template size<M>());
                        return data.template get<types::vector1d_t<T>>();
                    }
                }
            }

            template<size_t M>
            static auto check_size(const size_t& n, const size_t& m) {
                utils::dynamic_assert(n == m, "Wrong dimension size at level ", M, ". Expected ", m, ". Got ", n);
            }

            template<size_t M>
            static auto check_array_size(const json& data, const size_t& m) {
                utils::dynamic_assert(data.is_array(), "Expected array at level ", M);
                check_size<M>(data.size(), m);
            }

            template<size_t M>
            static auto read_data(const utils::dimensions<D...>& dims, const std::filesystem::path& path, const std::string& filename, const bool& binary) {
                const auto file_path = utils::make_file_path(path, filename);
                if (binary)
                    return ample::vector_reader<T>::template binary_read<M, D...>(std::ifstream(file_path, std::ios::binary), dims);
                return ample::vector_reader<T>::template read<M, D...>(std::ifstream(file_path), dims);
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

    }// namespace _impl;

#define CONFIG_DATA_FIELD(field, type)                                      \
    private:                                                                \
        mutable std::optional<type> _data_##field;                          \
    public:                                                                 \
        const type& field () const {                                        \
            if (!_data_##field.has_value()) {                               \
            	utils::dynamic_assert(_data.contains(#field),               \
                    "Missing field " #field);                               \
                _data_##field = type(_data[#field].template get<type>());   \
            }                                                               \
            return _data_##field.value();                                   \
        }                                                                   \
        bool has_##field() const {                                          \
            return _data_##field.has_value();                               \
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

#define CONFIG_INPUT_DATA(field, op, ...)                                   \
    private:                                                                \
        std::optional<__VA_ARGS__> _##field;                                \
    public:                                                                 \
        const auto& field () const {                                        \
            utils::dynamic_assert(_##field.has_value(),                     \
                "Field " #field " has no value");                           \
            return _##field.value() op;                                     \
        }                                                                   \
        bool has_##field () const {                                         \
            return _##field.has_value();                                    \
        }                                                                   \
        template<typename... Args>                                          \
        std::enable_if_t<0 < sizeof...(Args)> field(Args&&... args) {       \
            _##field = __VA_ARGS__(std::forward<Args>(args)...);            \
        }

#define CONFIG_INPUT_MULTI_DATA_DEF(field, ...)                             \
    private:                                                                \
        types::multi_optional<__VA_ARGS__> _##field;

#define CONFIG_INPUT_MULTI_DATA(name, field, op, ...)                       \
    public:                                                                 \
        const auto& name () const {                                         \
            if (!_##field.template has_value<__VA_ARGS__>())                \
                throw std::runtime_error("Field " #field " has no value");  \
            return _##field.template value<__VA_ARGS__>() op;               \
        }                                                                   \
        bool has_##name () const {                                          \
            return _##field.template has_value<__VA_ARGS__>();              \
        }                                                                   \
        template<typename... Args>                                          \
        std::enable_if_t<0 < sizeof...(Args)> name(Args&&... args) {        \
            _##field = __VA_ARGS__(std::forward<Args>(args)...);            \
        }

#define ASSERT_NO_VALUE(name, var, ...)                                     \
    utils::dynamic_assert(!var.__VA_ARGS__(),                               \
        "Multiple values provided for " name);

#define READ_INPUT_DATA(type, name, var, func, data, path, ...)             \
    if (type == name) {                                                     \
        __VA_ARGS__;                                                        \
        var = func(data, path);                                             \
        continue;                                                           \
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

            auto cb =
                    [keys=types::vector1d_t<std::unordered_set<std::string>>(1), path=std::vector<std::string>(1), indx=std::vector<int>(1, -2)]
                    (const int& depth, const json::parse_event_t event, json& parsed) mutable {
                        switch (event) {
                            case json::parse_event_t::array_end:
                                indx[depth + 1] = -2;
                            case json::parse_event_t::object_end:
                                keys[depth + 1].clear();
                                path[depth + 1].clear();
                                break;
                            case json::parse_event_t::array_start:
                            case json::parse_event_t::object_start:
                                if (keys.size() == depth + 1) {
                                    keys.emplace_back();
                                    path.emplace_back();
                                    indx.emplace_back(-2);
                                }

                                if (event == json::parse_event_t::array_start)
                                    indx[depth + 1] = -1;
                            case json::parse_event_t::value:
                                if (indx[depth] >= -1)
                                    ++indx[depth];
                                break;
                            case json::parse_event_t::key:
                            {
                                const auto p = keys[depth].insert(path[depth] = parsed.get<std::string>());
                                if (!p.second) {
                                    for (int i = 0; i < depth; ++i)
                                        if (indx[i] >= 0)
                                            path[i] = std::to_string(indx[i]);

                                    throw std::runtime_error(
                                            utils::join("Duplicate key \"",
                                                utils::join_it(path.begin() + 1, path.begin() + depth + 1, "/"),'"'));
                                }
                                break;
                            }
                            default:
                                break;
                        }

                        return true;
            };
            const json data = json::parse(std::ifstream(filename), cb, true, true);
            _data.merge_patch(data);

            const auto bls = bottom_layers().size();
            utils::dynamic_assert(betas().size() == n_layers() + bls,
                                  "Size of \"betas\"(", betas().size(), ") must be equal to n_layers + size of \"bottom_layers\"(", n_layers() + bls, ")");
            utils::dynamic_assert(bls == bottom_rhos().size() && bls == bottom_c1s().size() && bls == bottom_c2s().size(),
                                  "Size of \"bottom_layers\"(", bls, "), \"bottom_rhos\"(", bottom_rhos().size(),
                                  "), \"bottom_c1s\"(", bottom_c1s().size(), ") and \"bottom_c2s\"(", bottom_c2s().size(), ") must be the same");

            for (const auto& it : _data["input_data"]) {
                const auto type = it["type"].template get<std::string>();

                try {
                    READ_INPUT_DATA(type, "k0",          _k0,          _read_k1d_data,   it, _path, ASSERT_NO_VALUE("k0",          _k0,          has_value))
                    READ_INPUT_DATA(type, "phi_s",       _phi_s,       _read_k1d_data,   it, _path, ASSERT_NO_VALUE("phi_s",       _phi_s,       has_value))
                    READ_INPUT_DATA(type, "phi_j",       _phi_j,       _read_phi_j,      it, _path, ASSERT_NO_VALUE("phi_j",       _phi_j,       has_value))
                    READ_INPUT_DATA(type, "frequencies", _frequencies, _read_1d_data,    it, _path, ASSERT_NO_VALUE("frequencies", _frequencies, has_value))
                    READ_INPUT_DATA(type, "bathymetry",  _bathymetry,  _read_bathymetry, it, _path, ASSERT_NO_VALUE("bathymetry",  _bathymetry,  has_value))
                    READ_INPUT_DATA(type, "hydrology",   _hydrology,   _read_hydrology,  it, _path, ASSERT_NO_VALUE("hydrology",   _hydrology,   has_value))
                    READ_INPUT_DATA(type, "receivers",   _receivers,   _read_receivers,  it, _path, ASSERT_NO_VALUE("receivers",   _receivers,   has_value))

                    READ_INPUT_DATA(type, "k_j", _k_j, _read_k_j<T>, it, _path,
                                      ASSERT_NO_VALUE("k_j", _k_j, template has_value<types::vector1d_t<utils::linear_interpolated_data_2d<T, T>>>))
                    READ_INPUT_DATA(type, "complex_k_j", _k_j, _read_k_j<std::complex<T>>, it, _path,
                                      ASSERT_NO_VALUE("complex k_j", _k_j, template has_value<types::vector1d_t<utils::linear_interpolated_data_2d<T, T>>>))

                    READ_INPUT_DATA(type, "source_function", (std::tie(_times, _source_function)), _read_source_function<T>, it, _path,
                                      ASSERT_NO_VALUE("times", _times, has_value); ASSERT_NO_VALUE("source function", _source_function, has_value))
                    READ_INPUT_DATA(type, "source_spectrum", (std::tie(_frequencies, _source_spectrum)), _read_source_function<std::complex<T>>, it, _path,
                                      ASSERT_NO_VALUE("frequencies", _frequencies, has_value); ASSERT_NO_VALUE("source spectrum", _source_spectrum, has_value))

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
        CONFIG_DATA_FIELD(z0, T)
        CONFIG_DATA_FIELD(z1, T)
        CONFIG_DATA_FIELD(nz, size_t)
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

        CONFIG_INPUT_DATA(bathymetry, [0], utils::linear_interpolated_data_2d<T>)
        CONFIG_INPUT_DATA(hydrology, [0], utils::delaunay_interpolated_data_2d<T>)
        CONFIG_INPUT_DATA(receivers, ;, types::vector1d_t<types::point<T>>)
        CONFIG_INPUT_DATA(source_function, ;, types::vector1d_t<T>)
        CONFIG_INPUT_DATA(source_spectrum, ;, types::vector1d_t<std::complex<T>>)
        CONFIG_INPUT_DATA(frequencies, ;, types::vector1d_t<T>)
        CONFIG_INPUT_DATA(times, ;, types::vector1d_t<T>)
        CONFIG_INPUT_DATA(k0, [_index], types::vector2d_t<T>)
        CONFIG_INPUT_DATA(phi_s, [_index], types::vector2d_t<T>)
        CONFIG_INPUT_DATA(phi_j, [_index], types::vector1d_t<utils::linear_interpolated_data_3d<T>>)

        CONFIG_INPUT_MULTI_DATA_DEF(k_j, types::vector1d_t<utils::linear_interpolated_data_2d<T>>, types::vector1d_t<utils::linear_interpolated_data_2d<T, std::complex<T>>>)
        CONFIG_INPUT_MULTI_DATA(k_j, k_j, [_index], types::vector1d_t<utils::linear_interpolated_data_2d<T>>)
        CONFIG_INPUT_MULTI_DATA(complex_k_j, k_j, [_index], types::vector1d_t<utils::linear_interpolated_data_2d<T, std::complex<T>>>)

        auto all_k0() const {
            return _k0.value();
        }

        auto all_phi_s() const {
            return _phi_s.value();
        }

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

        auto mnx() const {
            return _data.contains("mnx") ? _data["mnx"].get<size_t>() : bathymetry().x().size();
        }

        auto mnx(const size_t& n) {
            _data["mnx"] = n;
        }

        auto mny() const {
            return _data.contains("mny") ? _data["mny"].get<size_t>() : bathymetry().y().size();
        }

        auto mny(const size_t& n) {
            _data["mny"] = n;
        }

        auto mnz() const {
            return _data.contains("mnz") ? _data["mnz"].get<size_t>() : nz();
        }

        auto mnz(const size_t& n) {
            _data["mnz"] = n;
        }

        template<typename V = T>
        auto create_modes(const size_t& c = 0, const bool show_progress = false) const {
            if (_k_j.template has_value<types::vector1d_t<utils::linear_interpolated_data_2d<T, V>>>() &&
                _phi_j.has_value()) {
                const auto& k_j = std::get<types::vector1d_t<utils::linear_interpolated_data_2d<T, V>>>(_k_j)[_index];
                const auto& phi_j = _phi_j.value()[_index];
                utils::dynamic_assert(k_j.size() >= c,
                      "Insufficient number of modes. Expect no less than ", c, ", but got ", k_j.size());
                utils::dynamic_assert(phi_j.size() >= c,
                      "Insufficient number of modes. Expect no less than ", c, ", but got ", phi_j.size());
                return std::make_tuple(k_j, phi_j);
            }

            const auto xn = mnx();
            const auto yn = mny();

            modes<T, V> modes(*this, utils::mesh_1d(z0(), z1(), mnz()));
            return modes.interpolated_field(xn, yn, utils::progress_bar_callback(xn * yn, "Modes", show_progress), c);
        }

        template<typename V = T>
        auto create_const_modes(const size_t& c = 0, const bool show_progress = false) const {
            if (_k_j.template has_value<types::vector1d_t<utils::linear_interpolated_data_2d<T, V>>>() &&
                _phi_j.has_value()) {
                const auto& k_j = std::get<types::vector1d_t<utils::linear_interpolated_data_2d<T, V>>>(_k_j)[_index];
                const auto& phi_j = _phi_j.value()[_index];
                utils::dynamic_assert(k_j.size() >= c,
                                      "Insufficient number of wave numbers. Expected no less than ", c, ", but got ", k_j.size());
                utils::dynamic_assert(phi_j.size() >= c,
                                      "Insufficient number of modal functions. Expected no less than ", c, ", but got ", phi_j.size());

                const auto size = std::min(k_j.size(), phi_j.size());

                types::vector2d_t<V> const_k_j(size);
                types::vector3d_t<T> const_phi_j(size);

                for (size_t j = 0; j < size; ++j) {
                    const_k_j[j] = k_j[j].data()[0];
                    const_phi_j[j] = phi_j[j].data()[0];
                }

                return std::make_tuple(
                    utils::linear_interpolated_data_1d<T, V>(k_j.template get<1>(), std::move(const_k_j)),
                    utils::linear_interpolated_data_2d<T, T>(phi_j.template get<1>(), phi_j.template get<2>(), std::move(const_phi_j))
                );
            }

            return create_const_modes<V>(utils::mesh_1d(z0(), z1(), mnz()), c, show_progress);
        }

        template<typename V = T>
        auto create_const_modes(const types::vector1d_t<T>& z, const size_t& c = 0, const bool show_progress = false) const {
            const auto yn = mny();

            modes<T, V> modes(*this, z);
            return modes.interpolated_line(x0(), yn, utils::progress_bar_callback(yn, "Modes", show_progress), c);
        }

        auto create_source_modes(const size_t& c = 0) const {
            if (has_k0() && has_phi_s())
                return std::make_tuple(k0(), phi_s());

            modes<T> modes(*this, { z_s() });
            const auto [k0, phi_s] = modes.point(0, y_s(), c);
            return std::make_tuple(k0, utils::make_vector(phi_s, [](const auto& data) { return data[0]; }));
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
                { "x0", T(0) },
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

                std::filesystem::copy(utils::make_file_path(path, std::filesystem::path(data.template get<std::string>())),
                    filename, std::filesystem::copy_options::overwrite_existing);

                data = filename.generic_string();
                return;
            }

            if (data.is_array()) {
                for (auto& it : data)
                    _copy_files(it, count, depth - 1, path, output, binary);

                return;
            }

            if (data.is_object())
                return;

            throw std::logic_error("Data must be either string, array or object");
        }

        static auto _read_hydrology(const json& data, const std::filesystem::path& path) {
            const auto [dimensions, inp_data] = _impl::input_data<T, T, T>(data, path);

            const auto& x = dimensions.template get<1>();
            const auto& z = dimensions.template get<0>();

            types::vector1d_t<T> px, pz, vv;
            px.reserve(z.size() * x.size());
            pz.reserve(z.size() * x.size());
            vv.reserve(z.size() * x.size());
            for (size_t i = 0; i < z.size(); ++i)
                for (size_t j = 0; j < x.size(); ++j)
                    if (inp_data[i][j] > -T(1e-10)) {
                        pz.emplace_back(z[i]);
                        px.emplace_back(x[j]);
                        vv.emplace_back(inp_data[i][j]);
                    }
            pz.shrink_to_fit();
            px.shrink_to_fit();
            vv.shrink_to_fit();

            return utils::delaunay_interpolated_data_2d<T>({px, pz}, vv);
        }

        static auto _read_bathymetry(const json& data, const std::filesystem::path& path) {
            const auto [dimensions, inp_data] = _impl::input_data<T, T, T>(data, path);

            return utils::linear_interpolated_data_2d<T>(
                dimensions.template get<0>(),
                dimensions.template get<1>(),
                inp_data
            );
        }

        static auto _read_receivers(const json& data, const std::filesystem::path& path) {
            return _impl::input_data<types::point<T>, utils::no_values_dim>(data, path).data;
        }

        template<typename V = T>
        static auto _read_source_function(const json& data, const std::filesystem::path& path) {
            const auto [dimensions, inp_data] = _impl::input_data<V, T>(data, path);
            return std::make_tuple(dimensions.template get<0>(), inp_data);
        }

        template<typename V = T>
        static auto _read_k_j(const json& data, const std::filesystem::path& path) {
            const auto [dimensions, inp_data] = _impl::input_data<V, utils::var_dim<utils::no_values_dim>, T, T>(data, path);
            return utils::make_vector(inp_data, [&dimensions](const auto& data) {
                return utils::linear_interpolated_data_2d<T, V>(
                    dimensions.template get<1>(),
                    dimensions.template get<2>(),
                    data
                );
            });
        }

        static auto _read_phi_j(const json& data, const std::filesystem::path& path) {
            const auto [dimensions, inp_data] = _impl::input_data<T, utils::var_dim<utils::no_values_dim>, T, T, T>(data, path);
            return utils::make_vector(inp_data, [&dimensions](const auto& data) {
                return utils::linear_interpolated_data_3d<T, T>(
                    dimensions.template get<1>(),
                    dimensions.template get<2>(),
                    dimensions.template get<3>(),
                    data
                );
            });
        }

        static auto _read_1d_data(const json& data, const std::filesystem::path& path) {
            return _impl::input_data<T, utils::no_values_dim>(data, path).data;
        }

        static auto _read_k1d_data(const json& data, const std::filesystem::path& path) {
            return _impl::input_data<T, utils::var_dim<utils::no_values_dim>>(data, path).data;
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

            throw std::runtime_error("Unknown coefficients type: " + type);
        }

    };

}// namespace ample
