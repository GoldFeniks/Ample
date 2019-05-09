#pragma once
#include <tuple>
#include <vector>
#include <istream>
#include "../utils/types.hpp"

namespace acstc {

    template<typename T = double>
    class table_reader {

    public:

        static auto read(std::istream& stream) {
            read_data data;
            std::tie(std::ignore, data.cols) = _read_line(stream);
            while (_find_something(stream)) {
                const auto [val, row] = _read_line(stream);
                if (row.size()) {
                    data.rows.push_back(std::move(val));
                    data.data.push_back(std::move(row));
                }
            }
            return data;
        }

    private:

        struct read_data {

            types::vector1d_t<T> cols, rows;
            types::vector2d_t<T> data;

        };

        static auto _read_line(std::istream& stream) {
            types::vector1d_t<T> data;
            T row, buff;
            stream >> row;
            while (_find_something(stream)) {
                stream >> buff;
                data.push_back(buff);
            }
            stream.get();
            return std::make_tuple(row, data);
        }

        static auto _find_something(std::istream& stream) {
            auto c = stream.peek();
            while (c != '\n' && c != EOF) {
                if (c == ' ' || c == '\t')
                    while (c == ' ' || c == '\t') {
                        stream.get();
                        c = stream.peek();
                    }
                else
                    return true;
            }
            return false;
        }

    };

}// namespace acstc
