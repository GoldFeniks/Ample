#pragma once
#include <tuple>
#include <string>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <filesystem>
#include "types.hpp"

namespace acstc {

    namespace utils {

        namespace __impl {

            template<size_t N, size_t M = N>
            struct gather_values {

                gather_values() = delete;

                template<typename Tuple, typename... Values>
                static auto get(const size_t& index, const Tuple& tuple, const Values&... values) {
                    return gather_values<N, M - 1>::get(
                            index, tuple, std::forward<const Values>(values)..., std::get<N - M>(tuple)[index]);
                }

            };

            template<size_t N>
            struct gather_values<N, 0> {

                gather_values() = delete;

                template<typename Tuple, typename... Values>
                static auto get(const size_t& index, const Tuple& tuple, const Values&... values) {
                    return std::make_tuple(std::forward<const Values>(values)...);
                }

            };

            template<typename It>
            class stride_iterator {

            public:

                explicit stride_iterator(It begin, It end, const size_t& k) : _it(begin), _end(end), _k(k) {}

                bool operator==(const stride_iterator& other) const {
                    return _it == other._it;
                }

                bool operator!=(const stride_iterator& other) const {
                    return !(*this == other);
                }

                bool operator==(const It& other) const {
                    return _it == other;
                }

                bool operator!=(const It& other) const {
                    return !(*this == other);
                }

                const auto& operator*() const {
                    return *_it;
                }

                stride_iterator& operator++() {
                    advance();
                    return *this;
                }

                const stride_iterator operator++(int) {
                    auto res = stride_iterator(_it, _end, _k);
                    advance();
                    return res;
                }

             private:

                void advance() {
                    if (std::distance(_it, _end) < _k)
                        _it = _end;
                    else
                        _it += _k;
                }

                It _it, _end;
                const size_t _k;

            };

        }// namespace __impl

        template<typename... T>
        class zip {

        public:

            explicit zip(const T&... values)
                : _values(std::forward<const T>(values)...), _size(std::get<0>(_values).size()) {}

            auto begin() const {
                return iterator(*this);
            }

            auto end() const {
                return iterator(*this, size());
            }

            auto size() const {
                return _size;
            }

            auto operator[](const size_t& index) const {
                return __impl::gather_values<N>::get(index, _values);
            }

        private:

            const std::tuple<const T&...> _values;
            const size_t _size;

            static constexpr auto N = sizeof...(T);

            class iterator {

            public:

                explicit iterator(const zip& owner, const size_t index = 0) : _owner(owner), _index(index) {}

                bool operator==(const iterator& other) const {
                    return _index == other._index;
                }

                bool operator!=(const iterator& other) const {
                    return !(*this == other);
                }

                auto operator*() const {
                    return _owner[_index];
                }

                iterator& operator++() {
                    ++_index;
                    return *this;
                }

                iterator operator++(int) {
                    auto res = iterator(_owner, _index);
                    ++_index;
                    return res;
                }

            private:

                const zip& _owner;
                size_t _index = 0;

            };

        };

        template<typename V, typename T>
        auto find_indices(const V& values, const T& value) {
            auto it = std::lower_bound(values.begin(), values.end(), value);
            if (it == values.end())
                return std::tuple<size_t, size_t>(0, 1);
            if (it == values.begin())
                ++it;
            const auto d = std::distance(values.begin(), it);
            return std::tuple<size_t, size_t>(d - 1, d);
        }

        template<typename T>
        auto mesh_1d(const T& a, const T& b, const size_t& n) {
            const auto h = (b - a) / (n - 1);
            types::vector1d_t<T> result;
            result.reserve(n);
            for (size_t i = 0; i < n; ++i)
                result.emplace_back(a + i * h);
            return result;
        }

        std::filesystem::path make_file_path(const std::filesystem::path root, const std::filesystem::path& filename) {
            if (!root.has_filename() || filename.is_absolute())
                return filename;
            std::filesystem::path result = root;
            return result.replace_filename(filename);
        }

        template<typename It>
        auto stride(It begin, It end, const size_t& k) {
            return std::make_tuple(
                __impl::stride_iterator(begin, end, k),
                __impl::stride_iterator(end,   end, k)
            );
        }

    }// namespace utils

}// namespace acstc