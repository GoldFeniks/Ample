#pragma once
#include <string>
#include <vector>
#include <cstddef>
#include <fstream>
#include <utility>
#include <type_traits>
#include <functional>
#include "types.hpp"

namespace acstc {

    namespace utils {

        namespace writer_bases {

            template<typename T>
            class writer_base {

            public:

                using value_type = T;

                virtual void before_write() {}
                virtual void after_write() {}
                virtual void write_one(const T&) = 0;

            };

            template<typename T>
            class stream_writer_base : public writer_base<T> {

            public:

                explicit stream_writer_base(const std::string& filename, std::string separator = " ", std::string ending = "\n",
                        std::ios_base::openmode mode = std::ios_base::out) :
                    _stream(filename, mode), _separator(std::move(separator)), _ending(std::move(ending)) {}

                void after_write() override {
                    _stream << _ending;
                }

                void write_one(const T& value) override {
                    _stream << value << _separator;
                }

                auto stream() {
                    return _stream;
                }

            private:

                std::ofstream _stream;
                const std::string _separator;
                const std::string _ending;

            };

            template<typename T>
            class binary_writer_base : public writer_base<T> {

            public:

                explicit binary_writer_base(const std::string& filename) : _stream(filename, std::ios_base::binary) {}

                void write_one(const T& value) override {
                    _stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
                }

                auto stream() {
                    return _stream;
                }

            private:

                std::ofstream _stream;

            };

            template<typename T, typename Base>
            class converter_writer_base : public writer_base<T> {

            public:

                template<typename F, typename... Args>
                explicit converter_writer_base(const F& convert, Args&&... args) :
                    _convert([convert](const T& value) { return convert(value); }), _base(std::forward<Args>(args)...) {}

                void before_write() {
                    _base.before_write();
                }

                void after_write() {
                    _base.after_write();
                }

                void write_one(const T& value) {
                    _base.write_one(_convert(value));
                }

            private:

                const std::function<typename Base::value_type(const T&)> _convert;
                Base _base;

            };

        }// namespace writer_bases

        template<typename T, typename Base>
        class writer : protected Base {

        public:

            using value_type = T;

            template<typename... Args>
            explicit writer(Args&&... args) : Base(std::forward<Args>(args)...) {}

            template<typename It>
            void write(It begin, It end) {
                this->before_write();
                while (begin != end)
                    this->write_one(*(begin++));
                this->after_write();
            }

            void write(const T* data, const size_t count) {
                write(data, data + count);
            }

            void write(const T& value) {
                this->before_write();
                this->write_one(value);
                this->after_write();
            }

            template<typename V, typename = std::enable_if_t<!std::is_same_v<T, V>>>
            void write(const V& data) {
                write(data.begin(), data.end());
            }

            void operator()(const T& value) {
                write(value);
            }

            template<typename V, typename = std::enable_if_t<!std::is_same_v<T, V>>>
            void operator()(const V& data) {
                write(data);
            }

        };

        template<typename T>
        class stream_writer : public writer<T, writer_bases::stream_writer_base<T>> {

        public:

            explicit stream_writer(const std::string& filename, const std::string& separator = " ",
                    const std::string& ending = "\n") :
                writer<T, writer_bases::stream_writer_base<T>>(filename, separator, ending) {}

        };

        template<typename T>
        class binary_writer : public writer<T, writer_bases::binary_writer_base<T>> {

        public:

            explicit binary_writer(const std::string& filename) :
                writer<T, writer_bases::binary_writer_base<T>>(filename) {}

        };

        template<typename T, typename Base>
        class converter_writer : public writer<T, writer_bases::converter_writer_base<T, Base>> {

        public:

            template<typename F, typename... Args>
            explicit converter_writer(const F& convert, Args&&... args) :
                    writer<T, writer_bases::converter_writer_base<T, Base>>(convert, std::forward<Args>(args)...) {}

        };

    };// namespace utils

}// namespace acstc
