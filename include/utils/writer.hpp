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

        template<typename T>
        struct basic_writer {

            virtual bool write(const T*, size_t) = 0;
            virtual bool write(const types::vector1d_t<T>&) = 0;

        };

        template<typename T, typename Base>
        class writer : public basic_writer<T>, protected Base {

        public:

            template<typename... Args>
            explicit writer(const size_t step, Args&&... args) : _step(step), Base(std::forward<Args>(args)...) {}

            template<typename It>
            bool write(It begin, It end) {
                if (_num_written++ % _step)
                    return false;
                this->before_write();
                while (begin != end)
                    this->write_one(*(begin++));
                this->after_write();
                return true;
            }

            bool write(const T* data, const size_t count) override {
                return write(data, data + count);
            };

            bool write(const types::vector1d_t<T>& data) override {
                return write(data.cbegin(), data.cend());
            }

        protected:

            const size_t _step;
            size_t _num_written = 0;

        };

        template<typename T>
        class text_writer : public writer<T, writer_bases::text_writer_base<T>> {

        public:

            explicit text_writer(const std::string& filename, const size_t step = 1, const std::string& separator = " ") :
                writer<T, writer_bases::text_writer_base<T>>(step, filename, separator) {}

        };

        template<typename T, typename Base>
        class converter_writer : public writer<T, writer_bases::converter_writer_base<T, Base>> {

        public:

            template<typename F, typename... Args>
            explicit converter_writer(const F& convert, const size_t step = 1, Args&&... args) :
                    writer<T, writer_bases::converter_writer_base<T, Base>>(step, convert, std::forward<Args>(args)...) {}

        };

    };// namespace utils

}// namespace acstc
