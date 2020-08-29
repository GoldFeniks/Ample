#pragma once
#include <array>
#include <complex>
#include <istream>
#include "nlohmann/json.hpp"

template<typename T, size_t N>
std::istream& operator>>(std::istream& stream, std::array<T, N>& array) {
    for (size_t i = 0; i < N; ++i)
        stream >> array[i];
    return stream;
}

template<typename T>
std::istream& operator>>(std::istream& stream, std::complex<T>& value) {
    T real, imag;
    stream >> real >> imag;
    value = std::complex<T>(real, imag);
    return stream;
}

namespace nlohmann {

    template<typename T>
    struct adl_serializer<std::complex<T>> {

        static void from_json(const nlohmann::json& data, std::complex<T>& value) {
            if (data.is_number()) {
                value.real(data.template get<T>());
                return;
            }

            if (data.is_array()) {
                if (data.size() != 2)
                    throw std::logic_error("Exactly two numbers must be provided for complex values");
                value = std::complex<T>(data[0].template get<T>(), data[1].template get<T>());
                return;
            }

            if (data.is_object()) {
                value = std::complex<T>(data.at("real").template get<T>(), data.at("imag").template get<T>());
                return;
            }

            throw std::logic_error("Cannot parse complex value");
        }

    };

}// namespace nlohmann
