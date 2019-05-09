#pragma once
#include <tuple>
#include <vector>
#include <istream>
#include "io/reader.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/interpolation.hpp"

namespace acstc {

    template<typename T = types::real_t, typename I = utils::linear_interpolation>
    class bathymetry {

    public:

        bathymetry() = delete;

        static auto from_binary(std::istream& stream) {
            uint32_t n, m;
            stream.read(reinterpret_cast<char*>(&n), sizeof(n));
            stream.read(reinterpret_cast<char*>(&m), sizeof(m));
            types::vector1d_t<T> x(n), y(m);
            types::vector2d_t<T> depths(n, types::vector1d_t<T>(m));
            stream.read(reinterpret_cast<char*>(x.data()), sizeof(T) * n);
            stream.read(reinterpret_cast<char*>(y.data()), sizeof(T) * m);
            for (auto& it : depths)
                stream.read(reinterpret_cast<char*>(it.data()), sizeof(T) * m);
            return utils::interpolated_data_2d<T, T, I>(std::move(x), std::move(y), std::move(depths));
        }

        static auto from_text(std::istream& stream) {
            auto data = table_reader<T>::read(stream);
            return utils::interpolated_data_2d<T, T, I>(std::move(data.rows), std::move(data.cols), std::move(data.data));
        }

    };

}// namespace acstc
