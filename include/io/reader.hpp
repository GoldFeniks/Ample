#pragma once
#include <tuple>
#include <vector>
#include <istream>
#include <sstream>
#include "convertors.hpp"
#include "../utils/types.hpp"
#include "../utils/utils.hpp"
#include "../utils/assert.hpp"
#include "../utils/dimensions.hpp"

namespace ample {

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
                                                    "Incorrect number of elements on line ", line, ". Expected ",
                                                    dims.template size<M>(i), ", but got ", result.back().size());
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
                                                "Incorrect number of elements on line ", line, ". Expected ",
                                                dims.template size<M>(), ", but got ", result.size());
                        ++line;
                    }

                    return result;
                }
        }

    }// namespace _impl

    template<typename T>
    struct reader {

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

}// namespace ample
