#pragma once
#include <tuple>
#include <utility>
#include <iostream>
#include "types.hpp"

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

            template<typename T, typename DCallback, typename CCallback>
            void ekc_callback(const T& data, const size_t k, size_t& ck, DCallback& data_callback, CCallback& count_callback) {
                if (ck % k == 0) {
                    data_callback(data);
                    count_callback(ck);
                }
                ++ck;
            }

        }// namespace __impl

        const auto& nothing_callback() {
            static auto callback = [](auto){};
            return callback;
        }

        //each kth call callback
        template<typename DCallback, typename CCallback>
        auto ekc_callback(const size_t k, DCallback& data_callback, CCallback& count_callback) {
            return [k, &data_callback, &count_callback, ck=size_t(0)](const auto& data) mutable {
                __impl::ekc_callback(data, k, ck, data_callback, count_callback);
            };
        }

        //each kth call callback
        template<typename DCallback, typename CCallback>
        auto ekc_callback(const size_t k, DCallback&& data_callback, CCallback&& count_callback) {
            return [k, ck=size_t(0), data_callback=std::forward<DCallback>(data_callback),
                    count_callback=std::forward<CCallback>(count_callback)]
                    (const auto& data) mutable {
                __impl::ekc_callback(data, k, ck, data_callback, count_callback);
            };
        }

        template<typename DCallback>
        auto ekc_callback(const size_t k, const DCallback& data_callback) {
            return ekc_callback(k, std::forward<DCallback>(data_callback), nothing_callback());
        }

        template<typename DCallback>
        auto ekc_callback(const size_t k, DCallback&& data_callback) {
            return ekc_callback(k, std::forward<DCallback>(data_callback), nothing_callback());
        }

        template<typename DCallback>
        auto progress_callback(const size_t k, DCallback& data_callback) {
            static auto callback = [](const size_t& n) { std::cout << n << std::endl; };
            return ekc_callback(k, std::forward<DCallback>(data_callback), callback);
        }

        template<typename DCallback>
        auto progress_callback(const size_t k, DCallback&& data_callback) {
            return ekc_callback(k, std::forward<DCallback>(data_callback),
                    [](const size_t& n) { std::cout << n << std::endl; });
        }

        auto progress_callback(const size_t k) {
            return progress_callback(k, nothing_callback());
        }

        template<typename... Callbacks>
        static auto callbacks(Callbacks&... callbacks) {
            return [callbacks=std::tie(std::forward<Callbacks>(callbacks)...)]
                    (const auto& data) mutable {
                __impl::caller<sizeof...(Callbacks) - 1>::call(callbacks, data);
            };
        }

        template<typename... Callbacks>
        static auto callbacks(Callbacks&&... callbacks) {
            return [callbacks=std::make_tuple(std::forward<Callbacks>(callbacks)...)]
                (const auto& data) mutable {
                __impl::caller<sizeof...(Callbacks) - 1>::call(callbacks, data);
            };
        }

        template<typename Callback>
        static auto callback_factory(Callback& callback) {
            return [&callback]() { return callbacks(callback); };
        }

        template<typename Callback>
        static auto callback_factory(Callback&& callback) {
            return [callback]() { return callbacks(callback); };
        }

        template<typename Callback>
        static auto callback_factory(const size_t k, Callback& callback) {
            return [k, &callback]() { return ekc_callback(k, callback); };
        }

        template<typename Callback>
        static auto callback_factory(const size_t k, Callback&& callback) {
            return [k, callback]() { return ekc_callback(k, callback); };
        }

        template<typename Callback>
        static auto progress_callback_factory(const size_t k, Callback& callback) {
            return [k, &callback]() { return progress_callback(k, callback); };
        }

        template<typename Callback>
        static auto progress_callback_factory(const size_t k, Callback&& callback) {
            return [k, callback]() { return progress_callback(k, callback); };
        }


    }// namespace utils

}// namespace acstc
