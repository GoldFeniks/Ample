#pragma once
#include <istream>
#include "io/reader.hpp"
#include "utils/types.hpp"
#include "utils/interpolation.hpp"

namespace acstc {

    template<typename T = types::real_t>
    class bathymetry {

    public:

        bathymetry() = delete;

        static auto from_table(types::vector1d_t<T>&& x, types::vector1d_t<T>&& y, types::vector2d_t<T>&& v) {
            return utils::linear_interpolated_data_2d<T>(x, y, v);
        }

        static auto from_table(const types::vector1d_t<T>& x, const types::vector1d_t<T>& y, const types::vector2d_t<T>& v) {
            return utils::linear_interpolated_data_2d<T>(x, y, v);
        }

        static auto from_binary(std::istream& stream) {
            auto [y, x, depths] = binary_table_reader<T>::read(stream);
            return from_table(std::move(x), std::move(y), std::move(depths));
        }

        static auto from_binary(std::ifstream&& stream) {
            return from_binary(stream);
        }

        static auto from_text(std::istream& stream) {
            const auto [y, x, depths] = table_reader<T>::read(stream);
            return utils::linear_interpolated_data_2d<T>(std::move(x), std::move(y), std::move(depths));
        }

        static auto from_text(std::istream&& stream) {
            return from_text(stream);
        }

    };

}// namespace acstc
