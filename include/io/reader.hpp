#pragma once
#include <tuple>
#include <vector>
#include <istream>
#include "../utils/types.hpp"
#include "../utils/utils.hpp"
#include "../utils/convertors.hpp"
#include "../utils/dimensions.hpp"

namespace acstc {

    namespace __impl {

        template<typename... T>
        auto find_something(std::istream& stream, const T&... chars) {
            auto c = stream.peek();
            while (((c != chars) && ...)) {
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

        template<typename T, typename V>
        auto read_line(std::istream& stream) {
            types::vector1d_t<V> data;
            T row;
            V buff;
            stream >> row;
            while (find_something(stream, '\r', '\n', EOF)) {
                stream >> buff;
                data.push_back(buff);
            }
            stream.get();
            return std::make_tuple(row, data);
        }


        template<typename T, size_t M, bool B, typename... D>
        auto read_vector(std::istream& stream, const utils::dimensions<D...>& dims) {
            if constexpr (M + 1 < sizeof...(D))
                if constexpr (utils::dimensions<D...>::template is_variable_dim<M>)
                    return utils::make_vector_i(dims.template size<M>(),
                        [&stream, &dims](const size_t& i) mutable {
                            return utils::make_vector_i(dims.template size<M>(i),
                                [&stream, &dims](const auto&) {
                                    return read_vector<T, M + 1, B, D...>(stream, dims);
                                }
                            );
                        }
                    );
                else
                    return utils::make_vector_i(dims.template size<M>(),
                        [&stream, &dims](const auto&) {
                            return read_vector<T, M + 1, B, D...>(stream, dims);
                        }
                    );
            else
                if constexpr (utils::dimensions<D...>::template is_variable_dim<M>) {
                    types::vector2d_t<T> result;

                    for (size_t i = 0; i < dims.template size<M>(); ++i) {
                        result.emplace_back(dims.template size<M>(i));
                        if constexpr (B)
                            stream.read(reinterpret_cast<char*>(result.back().data()), sizeof(T) * dims.template size<M>(i));
                        else
                            for (size_t j = 0; j < dims.template size<M>(i); ++j)
                                stream >> result.back()[j];
                    }

                    return result;
                } else {
                    types::vector1d_t<T> result(dims.template size<M>());

                    if constexpr (B)
                        stream.read(reinterpret_cast<char*>(result.data()), sizeof(T) * dims.template size<M>());
                    else
                        for (size_t j = 0; j < dims.template size<M>(); ++j)
                            stream >> result[j];

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
            std::tie(std::ignore, data.cols) = __impl::read_line<T, V>(stream);
            while (__impl::find_something(stream, EOF)) {
                auto [val, row] = __impl::read_line<T, V>(stream);
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

    template<typename T, typename V = T, typename S = uint32_t>
    class binary_table_reader {

    public:

        static auto read(std::istream& stream) {
            S n, m;
            stream.read(reinterpret_cast<char*>(&n), sizeof(S));
            stream.read(reinterpret_cast<char*>(&m), sizeof(S));
            read_data<T, V> data;
            data.data.resize(n, types::vector1d_t<V>(m));
            data.rows.resize(n);
            data.cols.resize(m);
            stream.read(reinterpret_cast<char*>(data.rows.data()), sizeof(T) * n);
            stream.read(reinterpret_cast<char*>(data.cols.data()), sizeof(T) * m);
            for (auto& it : data.data)
                stream.read(reinterpret_cast<char*>(it.data()), sizeof(V) * m);
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

            while (__impl::find_something(stream, EOF)) {
                stream >> a >> b;
                first.push_back(a);
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
            return __impl::read_vector<T, M, false, D...>(stream, dims);
        }

        template<size_t M = 0, typename... D>
        static auto binary_read(std::istream&& stream, const utils::dimensions<D...>& dims) {
            return __impl::read_vector<T, M, true, D...>(stream, dims);
        }

    };

}// namespace acstc
