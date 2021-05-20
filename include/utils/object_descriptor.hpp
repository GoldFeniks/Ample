#pragma once
#include <string>
#include <utility>
#include "assert.hpp"
#include "nlohmann/json.hpp"

namespace ample::utils {

    class object_descriptor {

    public:

        object_descriptor() = default;
        object_descriptor(std::string type, nlohmann::json parameters) : _type(std::move(type)), _parameters(std::move(parameters)) {}

        template<typename T, typename... Types, typename... Args>
        T construct(const Args&... args) const {
            static_assert(sizeof...(Types) == sizeof...(Args), "Incorrect number of arguments");
            (dynamic_assert(_parameters.contains(args), "Parameter \"", args, "\" is missing"), ...);
            return T(_parameters[args].template get<Types>()...);
        }

        template<typename... Args>
        bool has_args(const Args&... args) const {
            return (_parameters.contains(args) && ...);
        }

        [[nodiscard]] const auto& type() const {
            return _type;
        }

        [[nodiscard]] const auto& parameters() const {
            return _parameters;
        }

    private:

        std::string _type;
        nlohmann::json _parameters;

        friend struct nlohmann::adl_serializer<object_descriptor>;

    };

}// namespace ample::utils

namespace nlohmann {

    template<>
    struct adl_serializer<ample::utils::object_descriptor> {

        static void from_json(const json& data, ample::utils::object_descriptor& value) {
            ample::utils::dynamic_assert(data.contains("type"), "Missing \"type\" key from object description");
            ample::utils::dynamic_assert(data.contains("parameters"), "Missing \"parameters\" key from object description");
            value._type = data["type"].get<std::string>();
            value._parameters = data["parameters"];
        }

        static void to_json(json& json, const ample::utils::object_descriptor& value) {
            json = { { "type", value._type }, { "parameters", value._parameters } };
        }

    };

}// namespace nlohmann
