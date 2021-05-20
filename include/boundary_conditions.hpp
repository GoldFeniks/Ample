#pragma once
#include <cmath>
#include <tuple>
#include <string>
#include <cstddef>
#include <functional>
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/assert.hpp"
#include "nlohmann/json.hpp"
#include "utils/interpolation.hpp"

namespace ample {

    using nlohmann::json;

    template<typename T>
    struct pml_function {

        json description;
        std::function<T(const T&)> function;

        T operator()(const T& x) const {
            return function(x);
        }

    };

    template<typename A, typename V = std::complex<A>>
    class pml_boundary_conditions {

    public:

        using arg_t = A;
        using val_t = V;

        pml_boundary_conditions(const size_t& width, pml_function<A> function) : _width(width), _function(std::move(function)) {}

        template<typename VL>
        auto get_band_builder(const types::vector1d_t<VL>& k0, const types::vector2d_t<V>& b, const A& y0, const A& y1, const size_t& ny) const {
            return band_builder<VL>(*this, k0, b, y0, y1, ny);
        }

        auto width() const {
            return _width;
        }

    private:

        static constexpr auto im = V(0, 1);
        static constexpr auto ze = A(0);
        static constexpr auto on = A(1);
        static constexpr auto tw = A(2);

        size_t _width;
        pml_function<A> _function;

        template<typename VL>
        class band_builder {

        public:

            band_builder(const pml_boundary_conditions& owner, const types::vector1d_t<VL>& k0,
                         const types::vector2d_t<V>& b, const A& y0, const A& y1, const size_t& ny)
                         : _owner(owner), _nm(k0.size()), _ny(ny + 2 * owner._width), _nc(b[0].size()), _b(b), _k0(k0),
                           _ac(_nm, types::vector2d_t<V>(_nc, types::vector1d_t<V>(_ny))),
                           _bc(_nm, types::vector2d_t<V>(_nc, types::vector1d_t<V>(_ny))),
                           _cc(_nm, types::vector2d_t<V>(_nc, types::vector1d_t<V>(_ny))) {
                const auto hy = (y1 - y0) / (ny - 1);
                _y0 = y0 - _owner._width * hy;
                _y1 = y1 + _owner._width * hy;
                _sq_hy = std::pow(hy, 2);

                _c1.resize(owner._width);
                _c2.resize(owner._width);
                _c3.resize(owner._width);

                auto y = on;
                const auto dy = on / _owner._width;
                const auto d2 = dy / tw;
                for (size_t i = 0; i < owner._width; ++i, y -= dy) {
                    const auto aa = on + im * owner._function(y - d2);
                    const auto bb = on + im * owner._function(y);
                    const auto cc = on + im * owner._function(y + d2);

                    _c1[i] = on / (aa * bb);
                    _c2[i] = on / bb * (on / aa + on / cc) / _sq_hy;
                    _c3[i] = on / (cc * bb);
                }
            }

            void update(const types::vector2d_t<VL>& k) {
                update(k, 0, _nm);
            }

            void update(const types::vector2d_t<VL>& k, const size_t& j0, const size_t& j1) {
                for (size_t j = j0; j < j1; ++j) {
                    const auto sq_k0 = std::pow(_k0[j], 2);

                    for (size_t i = 0; i < _nc; ++i) {
                        const auto bk = _b[j][i] / sq_k0;
                        const auto dd = bk / _sq_hy;
                        const auto ty = tw / _sq_hy;

                        _ac[j][i].front() = _ac[j][i].back() = _cc[j][i].front() = _cc[j][i].back() = ze;
                        _bc[j][i].front() = _bc[j][i].back() = on;

                        size_t yi = 1;

                        const auto kf = std::pow(k[j].front(), 2);
                        for (size_t l = 1; l < _owner._width; ++l, ++yi) {
                            _ac[j][i][yi] = dd * _c1[l];
                            _bc[j][i][yi] = on + bk * (kf - sq_k0 - _c2[l]);
                            _cc[j][i][yi] = dd * _c3[l];
                        }

                        for (size_t l = 0; l < _ny - 2 * _owner._width; ++l, ++yi) {
                            _ac[j][i][yi] = _cc[j][i][yi] = dd;
                            _bc[j][i][yi] = on + bk * (std::pow(k[j][l], 2) - sq_k0 - ty);
                        }

                        const auto kb = std::pow(k[j].back(), 2);
                        for (size_t l = _owner._width - 1; l > size_t(0); --l, ++yi) {
                            _ac[j][i][yi] = dd * _c1[l];
                            _bc[j][i][yi] = on + bk * (kb - sq_k0 - _c2[l]);
                            _cc[j][i][yi] = dd * _c3[l];
                        }
                    }
                }
            }

            [[nodiscard]] auto coefficients() const {
                return std::tie(_ac, _bc, _cc);
            }

            [[nodiscard]] auto y0() const {
                return _y0;
            }

            [[nodiscard]] auto y1() const {
                return _y1;
            }

            [[nodiscard]] auto ny() const {
                return _ny;
            }

        private:

            const pml_boundary_conditions& _owner;

            A _sq_hy{}, _y0, _y1;
            const size_t _nm, _ny, _nc;
            const types::vector2d_t<V>& _b;
            const types::vector1d_t<VL>& _k0;

            types::vector1d_t<V> _c1, _c2, _c3;
            types::vector3d_t<V> _ac, _bc, _cc;

        };

    };

}// namespace ample

namespace nlohmann {

    template<typename T>
    struct adl_serializer<ample::pml_function<T>> {

        static void from_json(const nlohmann::json& data, ample::pml_function<T>& value) {
            const auto type = data["type"].get<std::string>();
            if (type == "cubic") {
                const auto a = data["/parameters/scale"_json_pointer].get<T>();
                value.function = [a](const auto& x) { return a * std::pow(x, 3); };
                value.description = { { "type", type }, { "parameters", { { "scale", a } } } };
                return;
            }

            if (type == "tabular") {
                const auto y = data["/parameters/values"_json_pointer].get<ample::types::vector1d_t<T>>();
                const auto x = ample::utils::mesh_1d<T>(0, 1, y.size());
                const auto interpolation = ample::utils::linear_interpolated_data_1d<T>(x, y);
                value.function = [interpolation](const auto& x) { return interpolation[0].point(x); };
                value.description = { { "type", type }, { "parameters", { { "values", y } } } };
                return;
            }

            ample::utils::dynamic_assert(false, "Unknown pml function type \"", type, '"');
        }

        static void to_json(nlohmann::json& json, const ample::pml_function<T>& value) {
            json = value.description;
        }

    };

}// namespace nlohmann
