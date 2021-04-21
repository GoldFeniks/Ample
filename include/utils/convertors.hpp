#pragma once
#include <array>
#include <regex>
#include <complex>
#include <istream>
#include <sstream>
#include "join.hpp"
#include "types.hpp"
#include "assert.hpp"
#include "nlohmann/json.hpp"

template<typename T>
std::istream& operator>>(std::istream& stream, ample::types::point<T>& point) {
    stream >> point.x >> point.y >> point.z;
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

            if (data.is_string()) {
                // Shamelessly stolen from https://stackoverflow.com/questions/50425322/c-regex-for-reading-complex-numbers
                constexpr auto raw_regex = R";(^(?=[iI.\d+-])([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?![iI.\d]))?([+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)?[iI])?$);";

                std::smatch match;
                const auto string = data.get<std::string>();
                std::regex regex(raw_regex, std::regex_constants::ECMAScript);

                ample::utils::dynamic_assert(std::regex_match(string, match, regex), "Couldn't match as complex value: ", string);

                const auto real = match[1].str();
                const auto imag = match[2].str();

                std::stringstream stream(
                        (real.empty() ? std::string("0") : real) + " " +
                        (imag.empty() ? std::string("0") : imag));
                stream >> value;
                return;
            }

            if (data.is_array()) {
                ample::utils::dynamic_assert(data.size() == 2, "Exactly two numbers must be provided for complex values");
                value = std::complex<T>(data[0].template get<T>(), data[1].template get<T>());
                return;
            }

            if (data.is_object()) {
                value = std::complex<T>(data.at("real").template get<T>(), data.at("imag").template get<T>());
                return;
            }

            throw std::runtime_error(ample::utils::join("Cannot parse complex value from ", data));
        }

    };

    template<typename T>
    struct adl_serializer<ample::types::point<T>> {

        static void from_json(const nlohmann::json& data, ample::types::point<T>& value) {
            if (data.is_object()) {
                ample::utils::dynamic_assert(
                    data.contains("x") && data.contains("y") && data.contains("z") && data.size() == 3 && data["x"].is_number() && data["y"].is_number() && data["z"].is_number(),
                    "Couldn't parse point value, expected an object in format { \"x\": <number>, \"y\": <number>, \"z\": <number> }, but got ", data
                );
                value = { data["x"].template get<T>(), data["y"].template get<T>(), data["z"].template get<T>() };
                return;
            }

            if (data.is_array()) {
                ample::utils::dynamic_assert(data.size() == 3, "Expected 3 values, but got ", data.size());
                value = { data[0].template get<T>(), data[1].template get<T>(), data[2].template get<T>() };
                return;
            }

            throw std::runtime_error(ample::utils::join("Cannot parse point value from ", data));
        }

        static void to_json(nlohmann::json& json, const ample::types::point<T>& value) {
            json = { { "x", value.x }, { "y", value.y }, { "z", value.z } };
        }

    };

}// namespace nlohmann
