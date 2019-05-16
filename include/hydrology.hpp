#pragma once
#include <istream>
#include <utility>
#include "io/reader.hpp"
#include "utils/types.hpp"
#include "utils/interpolation.hpp"

namespace acstc {

    template<typename T = types::real_t>
    class hydrology {

    public:

        hydrology() = delete;

        template<typename X, typename Z, typename V>
        static auto from_table(const X& x, const Z& z, const V& v) {
            types::vector1d_t<T> px, pz, vv;
            px.reserve(z.size() * x.size());
            pz.reserve(z.size() * x.size());
            vv.reserve(z.size() * x.size());
            for (size_t i = 0; i < x.size(); ++i)
                for (size_t j = 0; j < z.size(); ++j)
                    if (v[i][j] > -T(1e-10)) {
                        pz.emplace_back(z[j]);
                        px.emplace_back(x[i]);
                        vv.emplace_back(v[i][j]);
                    }
            pz.shrink_to_fit();
            px.shrink_to_fit();
            vv.shrink_to_fit();
            return utils::delaunay_interpolated_data_2d<T>({px, pz}, vv);
        }

        static auto from_text(std::istream& stream) {
            const auto [x, z, v] = table_reader<T>::read(stream);
            return from_table(x, z, v);
        }

        static auto from_text(std::istream&& stream) {
            return from_text(stream);
        }

        template<typename S = uint32_t>
        static auto from_binary(std::istream& stream) {
            S n, m;
            stream.read(reinterpret_cast<char*>(&n), sizeof(S));

            T x, z, v;
            types::vector1d_t<T> xs, zs, vs;
            for (size_t i = 0; i < n; ++i) {
                stream.read(reinterpret_cast<char*>(&x), sizeof(T));
                stream.read(reinterpret_cast<char*>(&m), sizeof(S));
                for (size_t j = 0; j < m; ++j) {
                    stream.read(reinterpret_cast<char*>(&z), sizeof(T));
                    stream.read(reinterpret_cast<char*>(&v), sizeof(T));
                    xs.push_back(x);
                    zs.push_back(z);
                    vs.push_back(v);
                }
            }
            xs.shrink_to_fit();
            zs.shrink_to_fit();
            vs.shrink_to_fit();
            return utils::delaunay_interpolated_data_2d<T>({xs, zs}, vs);
        }

        template<typename S = int32_t>
        static auto from_binary(std::istream&& stream) {
            return from_binary<S>(stream);
        }

    };

}// namespace acstc
