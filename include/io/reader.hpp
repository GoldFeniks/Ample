#pragma once
#include <tuple>
#include <vector>
#include <istream>
#include <sstream>
#include "../utils/types.hpp"
#include "../utils/utils.hpp"
#include "../utils/assert.hpp"
#include "../utils/convertors.hpp"
#include "../utils/dimensions.hpp"

namespace acstc {

    namespace _impl {

        template<typename T, typename V>
        auto read_table_line(std::istream& stream) {
            types::vector1d_t<V> data;

            std::string line;
            std::getline(stream, line);
            stream >> std::ws;

            if (line.empty())
                return std::make_tuple(T(), data);

            std::istringstream line_stream(line);

            T row;
            V buff;
            line_stream >> row;
            while (!line_stream.eof()) {
                line_stream >> buff;
                utils::dynamic_assert(!line_stream.fail(), "Incorrect data in file");

                data.push_back(buff);
                line_stream >> std::ws;
            }

            return std::make_tuple(row, data);
        }

        template<typename T>
        void read_line(std::istream& stream, types::vector1d_t<T>& data) {
            std::string line;
            std::getline(stream, line);
            stream >> std::ws;

            std::istringstream line_stream(line);
            T value;
            while (!line_stream.eof()) {
                line_stream >> value;
                utils::dynamic_assert(!line_stream.fail(), "Incorrect data in file");

                data.push_back(value);
                line_stream >> std::ws;
            }
        }


        template<typename T, size_t M, bool B, typename... D>
        auto read_vector(std::istream& stream, const utils::dimensions<D...>& dims, size_t& line) {
            if constexpr (M + 1 < sizeof...(D))
                if constexpr (utils::dimensions<D...>::template is_variable_dim<M>)
                    return utils::make_vector_i(dims.template size<M>(),
                        [&stream, &dims, &line](const size_t& i) mutable {
                            return utils::make_vector_i(dims.template size<M>(i),
                                [&stream, &dims, &line](const auto&) mutable {
                                    return read_vector<T, M + 1, B, D...>(stream, dims, line);
                                }
                            );
                        }
                    );
                else
                    return utils::make_vector_i(dims.template size<M>(),
                        [&stream, &dims, &line](const auto&) mutable {
                            return read_vector<T, M + 1, B, D...>(stream, dims, line);
                        }
                    );
            else
                if constexpr (utils::dimensions<D...>::template is_variable_dim<M>) {
                    types::vector2d_t<T> result;

                    for (size_t i = 0; i < dims.template size<M>(); ++i) {
                        result.emplace_back();
                        if constexpr (B) {
                            result.back().resize(dims.template size<M>(i));
                            stream.read(reinterpret_cast<char*>(result.back().data()), sizeof(T) * dims.template size<M>(i));
                            utils::dynamic_assert(stream.gcount() == sizeof(T) * dims.template size<M>(i), "Insufficient data in file");
                        }
                        else {
                            read_line(stream, result.back());
                            utils::dynamic_assert(result.back().size() == dims.template size<M>(i),
                                                  utils::join("Incorrect number of elements on line ", line, ". Expected ",
                                                              dims.template size<M>(i), ", but got ", result.back().size()));
                            ++line;
                        }
                    }

                    return result;
                } else {
                    types::vector1d_t<T> result;

                    if constexpr (B) {
                        result.resize(dims.template size<M>());
                        stream.read(reinterpret_cast<char*>(result.data()), sizeof(T) * dims.template size<M>());
                        utils::dynamic_assert(stream.gcount() == sizeof(T) * dims.template size<M>(), "Insufficient data in file");
                    }
                    else {
                        read_line(stream, result);
                        utils::dynamic_assert(result.size() == dims.template size<M>(),
                                              utils::join("Incorrect number of elements on line ", line, ". Expected ",
                                                          dims.template size<M>(), ", but got ", result.size()));
                        ++line;
                    }

                    return result;
                }
        }

    }// namespace __impl

    template<typename T, typename V = T>
    struct read_data {

        types::vector1d_t<T> cols, rows;
        types::vector2d_t<V> data;

    };

    template<typename T, typename V = T>
    class table_reader {

    public:

        static auto read(std::istream& stream) {
            read_data<T, V> data;
            std::tie(std::ignore, data.cols) = _impl::read_table_line<T, V>(stream);
            while (!(stream.eof() || stream.fail())) {
                auto [val, row] = _impl::read_table_line<T, V>(stream);
                if (row.size()) {
                    data.rows.push_back(std::move(val));
                    data.data.push_back(std::move(row));
                }
            }

            return data;
        }

        static auto read(std::istream&& stream) {
            return read(stream);
        }

    };

    template<typename T>
    void read_binary_values(std::istream& stream, T* data, const size_t& count) {
        stream.read(reinterpret_cast<char*>(data), sizeof(T) * count);
        utils::dynamic_assert(stream.gcount() == sizeof(T) * count, "Insufficient data in file");

    }

    template<typename T>
    void read_binary_values(std::istream&& stream, T* data, const size_t& count) {
        read_binary_values(stream, data, count);
    }

    template<typename T, typename V = T, typename S = uint32_t>
    class binary_table_reader {

    public:

        static auto read(std::istream& stream) {
            S n, m;
            read_binary_values(stream, &n, 1);
            read_binary_values(stream, &m, 1);

            read_data<T, V> data;
            data.data.resize(n, types::vector1d_t<V>(m));
            data.rows.resize(n);
            data.cols.resize(m);

            read_binary_values(stream, data.rows.data(), n);
            read_binary_values(stream, data.cols.data(), m);

            for (auto& it : data.data)
                read_binary_values(stream, it.data(), m);

            utils::dynamic_assert(stream.eof(), "Extra data in file");

            return data;
        }

        static auto read(std::istream&& stream) {
            return read(stream);
        }

    };

    template<typename T, typename V = T>
    class pairs_reader {

    public:

        static auto read(std::istream& stream) {
            types::vector1d_t<T> first;
            types::vector1d_t<V> second;
            T a;
            V b;

            while (!stream.eof()) {
                stream >> a >> std::ws;
                utils::dynamic_assert(!stream.fail(), "Incorrect data in file");
                utils::dynamic_assert(!stream.eof(), "Insufficient data in file");
                first.push_back(a);

                stream >> b;
                utils::dynamic_assert(!stream.fail(), "Incorrect data in file");
                second.push_back(b);
            }

            return std::make_tuple(first, second);
        }

        static auto read(std::istream&& stream) {
            return read(stream);
        }

    };

    template<typename T, typename V = T, typename S = uint32_t>
    class binary_pairs_reader {

    public:

        static auto read(std::istream& stream) {
            S n;
            stream.read(reinterpret_cast<char*>(&n), sizeof(S));

            types::vector1d_t<T> first(n);
            types::vector1d_t<V> second(n);

            for (size_t i = 0; i < n; ++i) {
                stream.read(reinterpret_cast<char*>(&first[i]), sizeof(T));
                stream.read(reinterpret_cast<char*>(&second[i]), sizeof(V));
            }

            return std::make_tuple(first, second);
        }

        static auto read(std::istream&& stream) {
            return read(stream);
        }

    };

    template<typename T>
    struct vector_reader {

        template<size_t M = 0, typename... D>
        static auto read(std::istream&& stream, const utils::dimensions<D...>& dims) {
            size_t line = 1;
            auto result = _impl::read_vector<T, M, false, D...>(stream, dims, line);
            stream >> std::ws;
            utils::dynamic_assert(stream.eof(), "Extra data in file");
            return result;
        }

        template<size_t M = 0, typename... D>
        static auto binary_read(std::istream&& stream, const utils::dimensions<D...>& dims) {
            size_t line = 0;
            auto result = _impl::read_vector<T, M, true, D...>(stream, dims, line);
            stream.get();
            utils::dynamic_assert(stream.eof(), "Extra data in file");
            return result;
        }

    };

}// namespace acstc
