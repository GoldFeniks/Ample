#pragma once
#include <tuple>
#include <string>
#include <cstddef>
#include <fstream>
#include "series.hpp"
#include "hydrology.hpp"
#include "bathymetry.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

namespace acstc {

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

        explicit config(const std::string& filename) : config() {
            std::ifstream in(filename);
            json data;
            in >> data;
            _data.merge_patch(data);
            _bathymetry = _create_bathymetry(_data["bathymetry"]);
            _hydrology = _create_hydrology(_data["hydrology"]);
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
        CONFIG_DATA_FIELD(zr, T)

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

    private:

        json _data;
        utils::linear_interpolated_data_2d<T> _bathymetry;
        utils::delaunay_interpolated_data_2d<T> _hydrology;
        T _a, _b, _c;

        static json _default_data() {
            return {
                { "mode_subset", -1 },
                { "ppm", size_t(2) },
                { "ordRich", size_t(3) },
                { "f", T(100) },
                { "zr", T(30) },
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
                { "x0", T(0) },
                { "x1", T(25000) },
                { "nx", size_t(2001) },
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

        auto _fill_coefficients(const json& data) {
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
