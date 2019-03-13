#pragma once
#include <string>
#include <vector>
#include <cstddef>
#include <fstream>
#include <utility>
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
            class text_writer_base : public writer_base<T> {

            public:

                text_writer_base(const std::string& filename, std::string separator = " ") :
                    _stream(filename), _separator(std::move(separator)) {}

                void after_write() override {
                    _stream << '\n';
                }

                void write_one(const T& value) override {
                    _stream << value << _separator;
                }

            private:

                std::ofstream _stream;
                const std::string _separator;

            };

            template<typename T, typename Base>
            class converter_writer_base : public writer_base<T> {

            public:

                template<typename F, typename... Args>
                converter_writer_base(const F& convert, Args&&... args) :
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
            };

            void write(const types::vector1d_t<T>& data) {
                write(data.cbegin(), data.cend());
            }

            void write(const T& value) {
                this->before_write();
                this->write_one(value);
                this->after_write();
            }

            void operator()(const types::vector1d_t<T>& data) {
                write(data);
            }

            void operator()(const T& value) {
                write(value);
            }

        };

        template<typename T>
        class text_writer : public writer<T, writer_bases::text_writer_base<T>> {

        public:

            explicit text_writer(const std::string& filename, const std::string& separator = " ") :
                writer<T, writer_bases::text_writer_base<T>>(filename, separator) {}

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
