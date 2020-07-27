#pragma once
#include <tuple>
#include <utility>
#include <iostream>
#include <functional>
#include "types.hpp"
#include "verbosity.hpp"
#include "progress_bar.hpp"

namespace acstc {

    namespace utils {

        namespace __impl {

            template<size_t N>
            struct caller {

                caller() = delete;

                template<typename Tuple, typename... Args>
                static void call(Tuple& tuple, Args&&... args) {
                    caller<N - 1>::call(tuple, std::forward<Args>(args)...);
                    std::get<N>(tuple)(std::forward<Args>(args)...);
                }

            };

            template<>
            struct caller<0> {

                caller() = delete;

                template<typename Tuple, typename... Args>
                static void call(Tuple& tuple, Args&&... args) {
                    std::get<0>(tuple)(std::forward<Args>(args)...);
                }

            };

            template<typename>
            struct _empty {};

            template<typename T>
            auto propagate(T& value, _empty<T&>) {
                return std::reference_wrapper<T>(value);
            }

            template<typename T>
            auto&& propagate(T& value, _empty<T>) {
                return std::move(value);
            }

            template<typename DCallback, typename CCallback>
            class ekc_callback {

            public:

                ekc_callback(const size_t& k, DCallback&& data_callback, CCallback&& count_callback) : 
                   _k(k), _data_callback(std::move(data_callback)), _count_callback(std::move(count_callback)) {}

                ekc_callback(ekc_callback&& other) : 
                    _ck(std::move(other._ck)), _k(other._k),
                    _data_callback(std::move(other._data_callback)), _count_callback(std::move(other._count_callback)) {}

                template<typename T, typename D>
                void operator()(const T& x, const D& data) {
                    if (_ck % _k == 0) {
                        _data_callback(x, data);
                        _count_callback(_ck);
                    }
                    ++_ck;
                }

                void reset() {
                    _ck = 0;
                }

            private:

                size_t _ck = 0;
                const size_t _k;
                DCallback _data_callback;
                CCallback _count_callback;

            };

            class progress_bar_callback {

            public:

                progress_bar_callback(const size_t& n, const std::string& title, const bool& enabled = true, 
                    const progress_bar::on_end& leave = progress_bar::on_end::global) :
                    _create_bar([=]() {
                        return new progress_bar(n, title, enabled, leave);
                    })
                {}

                ~progress_bar_callback() {
                    if (_bar)
                        delete _bar;
                }

                template<typename... T>
                void operator()(const T&...) {
                    if (!_bar)
                        _bar = _create_bar();
                    _bar->next();
                }

            private:

                progress_bar* _bar = nullptr;
                const std::function<progress_bar*()> _create_bar;


            };

        }// namespace __impl

        const auto& nothing_callback() {
            static auto callback = [](const auto&...){};
            return callback;
        }

        //each kth call callback
        template<typename DCallback, typename CCallback>
        auto ekc_callback(const size_t k, DCallback&& data_callback, CCallback&& count_callback) {
            return __impl::ekc_callback(k, 
                __impl::propagate(data_callback, __impl::_empty<DCallback>{}),
                __impl::propagate(count_callback, __impl::_empty<CCallback>{}));
        }

        template<typename DCallback>
        auto ekc_callback(const size_t k, DCallback&& data_callback) {
            return ekc_callback(k, std::forward<DCallback>(data_callback), nothing_callback());
        }

        template<typename DCallback>
        auto progress_callback(const size_t k, const verbosity& verbosity, DCallback&& data_callback, const size_t& level = 2) {
            return ekc_callback(k, std::forward<DCallback>(data_callback),
                    [&verbosity, level](const size_t& n) { verbosity.verboseln(level, n); });
        }

        template<typename DCallback>
        auto progress_callback(const size_t k, DCallback&& data_callback, const size_t& level = 2) {
            return progress_callback(k, verbosity::instance(), std::move(data_callback), level);
        }

        auto progress_callback(const size_t k) {
            return progress_callback(k, nothing_callback());
        }

        template<typename... Callbacks>
        auto callbacks(Callbacks&&... callbacks) {
            return [callbacks=std::tuple<Callbacks...>(std::forward<Callbacks>(callbacks)...)]
                (auto&&... values) mutable {
                __impl::caller<sizeof...(Callbacks) - 1>::call(callbacks, std::forward<decltype(values)>(values)...);
            };
        }

        auto progress_bar_callback(const size_t& n, const std::string& title, const bool& enabled = true, 
            const progress_bar::on_end& leave = progress_bar::on_end::global) {
            return __impl::progress_bar_callback(n, title, enabled, leave);
        }

    }// namespace utils

}// namespace acstc
