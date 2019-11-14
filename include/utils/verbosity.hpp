#pragma once
#include <cstddef>
#include <iostream>

namespace acstc {

    namespace utils {

        namespace __impl {

            template<typename T, typename... V>
            struct verbose {

                static void call(const T& value, const V&... values) {
                    std::cout << value;
                    verbose<V...>::call(values...);
                }

            };

            template<typename T>
            struct verbose<T> {

                static void call(const T& value) {
                    std::cout << value;
                }

            };


        }// namespace __impl

        class verbosity {

        public:

            size_t level = 0;

            static verbosity& instance() {
                static verbosity verbosity;
                return verbosity;
            }

            template<typename... T>
            void verbose(const size_t level, const T&... values) const {
                if (level <= this->level)
                    __impl::verbose<T...>::call(values...);
            }

            template<typename... T>
            void verboseln(const size_t level, const T&... values) const {
                if (level <= this->level) {
                    verbose(level, values...);
                    std::cout << std::endl;
                }
            }

        };

        template<typename... T>
        void verbose(const size_t level, const T&... values) {
            verbosity::instance().verbose<T...>(level, values...);
        }

        template<typename... T>
        void verboseln(const size_t level, const T&... values) {
            verbosity::instance().verboseln<T...>(level, values...);
        }

    }// namespace utils

}// namespace acstc