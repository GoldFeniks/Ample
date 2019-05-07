#include <tuple>
#include <string>
#include <cstddef>
#include <fstream>
#include "bathymetry.hpp"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

namespace acstc {

#define CONFIG_DATA_FIELD(field, type) const auto field () const { return _data[#field].template get<type>(); }

    template<typename T = double>
    class config {

    public:

        config() : _data(_default_data()), _bathymetry(_create_bathymetry(_data["bathymetry"])) {}

        explicit config(const std::string& filename) : config() {
            std::ifstream in(filename);
            json data;
            in >> data;
            _data.merge_patch(data);
            _bathymetry = _create_bathymetry(_data["bathymetry"]);
        }

        CONFIG_DATA_FIELD(mode_subset, int)
        CONFIG_DATA_FIELD(x0, T)
        CONFIG_DATA_FIELD(x1, T)
        CONFIG_DATA_FIELD(nx, size_t)
        CONFIG_DATA_FIELD(y0, T)
        CONFIG_DATA_FIELD(y1, T)
        CONFIG_DATA_FIELD(ny, size_t)
        CONFIG_DATA_FIELD(ppm, size_t)
        CONFIG_DATA_FIELD(ordRich, size_t)
        CONFIG_DATA_FIELD(f, T)

        const auto data() const {
            return _data;
        }

        const auto x_bounds() const {
            return std::make_tuple(x0(), x1());
        }

        const auto y_bounds() const {
            return std::make_tuple(y0(), y1());
        }

        const auto bounds() const {
            return std::tuple_cat(x_bounds(), y_bounds());
        }

        const auto bathymetry() const {
            return _bathymetry;
        }

    private:

        json _data;
        ::acstc::bathymetry<T> _bathymetry;

        static json _default_data() {
            return {
                { "mode_subset", -1 },
                { "ppm", size_t(2) },
                { "ordRich", size_t(3) },
                { "f", T(100) },
                { "bathymetry",
                  { "values",
                    {
                      { "x", { T(0), T(25000) } },
                      { "y", { T(-3500), T(3500) } },
                      { "depths",
                        {
                          { T(25), T(375) },
                          { T(25), T(375) }
                        }
                      }
                    }
                  }
                },
                { "x0", T(0) },
                { "x1", T(2000) },
                { "nx", size_t(2001) },
                { "y0", -T(1000) },
                { "y1",  T(1000) },
                { "ny",  size_t(2001) }
            };
        }

        static auto _create_bathymetry(const json& data) {
            const auto type = data[0].template get<std::string>();
            if (type == "values") {
                return ::acstc::bathymetry<T>(data["/1/x"_json_pointer].template get<types::vector1d_t<T>>(),
                                              data["/1/y"_json_pointer].template get<types::vector1d_t<T>>(),
                                              data["/1/depths"_json_pointer].template get<types::vector2d_t<T>>());
            }
            if (type == "text_file") {
                std::ifstream in(data[1].template get<std::string>());
                return ::acstc::bathymetry<T>::from_text(in);
            }
            if (type == "binary_file") {
                std::ifstream in(data[1].template get<std::string>(), std::ios::binary);
                return ::acstc::bathymetry<T>::from_binary(in);
            }
            throw std::logic_error("Unknown bathymetry type: " + type);
        }

    };

}// namespace acstc
