#pragma once
#include <tuple>
#include <vector>
#include <istream>
#include "io/reader.hpp"
#include "utils/types.hpp"
#include "utils/interpolation.hpp"

namespace acstc {

    template<typename T = types::real_t, typename V = T, typename I = utils::linear_interpolation>
    class bathymetry {

    public:

        static constexpr auto eps = T(1e-10);

        bathymetry() = delete;

        bathymetry(const types::vector1d_t<T>& x, const types::vector1d_t<T>& y, const types::vector2d_t<V>& depths) :
            _x(x), _y(y), _depths(depths) {}

        bathymetry(types::vector1d_t<T>&& x, types::vector1d_t<T>&& y, types::vector2d_t<V>&& depths) :
            _x(std::move(x)), _y(std::move(y)), _depths(std::move(depths)) {}

        static auto from_binary(std::istream& stream) {
            uint32_t n, m;
            stream.read(reinterpret_cast<char*>(&n), sizeof(n));
            stream.read(reinterpret_cast<char*>(&m), sizeof(m));
            types::vector1d_t<T> x(n), y(m);
            types::vector2d_t<T> depths(n, types::vector1d_t<T>(m));
            stream.read(reinterpret_cast<char*>(x.data()), sizeof(T) * n);
            stream.read(reinterpret_cast<char*>(y.data()), sizeof(T) * m);
            for (auto& it : depths)
                stream.read(reinterpret_cast<char*>(it.data()), sizeof(V) * m);
            return bathymetry(std::move(x), std::move(y), std::move(depths));
        }

        static auto from_text(std::istream& stream) {
            auto data = table_reader<T>::read(stream);
            return bathymetry(std::move(data.rows), std::move(data.cols), std::move(data.data));
        }

        void remove_insignificant_rows(const T& tolerance) {
            size_t j = 0, c = 0;
            for (size_t i = 2; i < _depths.size(); ++i) {
                for (size_t k = 0; k < _depths[i].size(); ++k)
                    if (std::abs(I::template line_point(_depths[j][k], _depths[i][k], _x[j], _x[i], _x[j + 1]) - _depths[j + 1][k]) > tolerance) {
                        if (c > 0) {
                            std::swap(_depths[i], _depths[j + 2]);
                            std::swap(_x[i], _x[j + 2]);
                        }
                        j += 1;
                        goto next;
                    }
                std::swap(_depths[i], _depths[j + 1]);
                std::swap(_x[i], _x[j + 1]);
                c += 1;
                next:
                continue;
            }
            _depths.resize(_depths.size() - c);
            _x.resize(_x.size() - c);
        }

        auto depth(const T& x, const T& y) const {
            const auto [ix, jx] = find_indices(_x, x);
            const auto [iy, jy] = find_indices(_y, y);
            return I::template field_point(_depths[ix][iy], _depths[ix][jy], _depths[jx][iy], _depths[jx][jy],
                    _x[ix], _x[jx], _y[iy], _y[jy], x, y);
        }

        auto depths_at(const T& x, const T& y0, const T& y1, const size_t n) const {
            const auto [ix, jx] = find_indices(_x, x);
            if (std::abs(_x[ix] - x) <= eps)
                return I::template line(y0, y1, n, _y, _depths[ix]);
            if (std::abs(_x[jx] - x) <= eps)
                return I::template line(y0, y1, n, _y, _depths[jx]);
            return I::template field(x, x, 1, y0, y1, n, _x, _y, _depths)[0];
        }

        auto depths_at(const T& x, const size_t n) const {
            return depths_at(x, _y.front(), _y.back(), n);
        }

        auto depths(const T& x0, const T& x1, const size_t nx, const T& y0, const T& y1, const size_t ny) const {
            return I::template field(x0, x1, nx, y0, y1, ny, _x, _y, _depths);
        }

        auto depths(const size_t nx, const size_t ny) const {
            return depths(_x.front(), _x.back(), nx, _y.front(), _y.back(), ny);
        }

        auto depths() const {
            return _depths;
        }

        auto x() const {
            return _x;
        }

        auto y() const {
            return _y;
        }

    private:

        types::vector1d_t<T> _x, _y;
        types::vector2d_t<V> _depths;

        template<typename Val>
        static auto find_indices(const Val& values, const T& value) {
            auto it = std::lower_bound(values.begin(), values.end(), value);
            if (it == values.end())
                return std::tuple<size_t, size_t>(0, 1);
            if (it == values.begin())
                ++it;
            const auto d = std::distance(values.begin(), it);
            return std::tuple<size_t, size_t>(d - 1, d);
        }

    };

}// namespace acstc
