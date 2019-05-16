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

        static auto from_binary(std::istream& stream) {
            const auto [y, x, depths] = binary_table_reader<T>::read(stream);
            return utils::linear_interpolated_data_2d<T>(std::move(x), std::move(y), std::move(depths));
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
