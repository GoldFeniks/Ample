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

        }// namespace __impl

        auto nothing_callback() {
            return [](auto){};
        }

        //each kth call callback
        template<typename DCallback, typename CCallback>
        auto ekc_callback(const size_t k, DCallback&& data_callback, CCallback&& count_callback) {
            return [k, ck=0, data_callback=std::forward<DCallback>(data_callback),
                    count_callback=std::forward<CCallback>(count_callback)]
                    (const auto& data) mutable {
                if (ck % k == 0) {
                    data_callback(data);
                    count_callback(ck);
                }
                ++ck;
            };
        }

        template<typename DCallback>
        auto ekc_callback(const size_t k, DCallback&& data_callback) {
            return ekc_callback(k, std::forward<DCallback>(data_callback), nothing_callback());
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
        static auto callbacks(Callbacks&&... callbacks) {
            return [callbacks=std::make_tuple(std::forward<Callbacks>(callbacks)...)]
                (const auto& data) mutable {
                __impl::caller<sizeof...(Callbacks) - 1>::call(callbacks, data);
            };
        }

    }// namespace utils

}// namespace acstc
