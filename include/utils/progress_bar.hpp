#pragma once
#include <chrono>
#include <string>
#include <cstddef>
#include <iostream>
#include <algorithm>
#include "boost/lexical_cast.hpp"

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/ioctl.h>
#endif

namespace acstc {

    namespace utils {

        class progress_bar {

        public:

            enum class on_end {
                clear, leave, global
            };

            inline static bool clear_on_end = false;

            progress_bar(const size_t& n, const std::string& description = "", 
                const bool& enabled = true, const on_end& end = on_end::global) : 
                _end(end),
                _n(n),
                _description(description),
                _enabled(enabled)
            {
#ifdef _WIN32
                _enable_vt100();
#endif
                write_progress();
            }

            progress_bar(const progress_bar&) = delete;
            progress_bar(progress_bar&&) = delete;

            progress_bar& operator=(const progress_bar&) = delete;
            progress_bar& operator=(progress_bar&) = delete;

            void operator()() {
                next();
            }

            ~progress_bar() {
                if (!_first && (_end == on_end::clear || _end == on_end::global && clear_on_end))
                    clear();
                else
                    printf("\n");
            }

            void next() {
                if (_cur >= _n)
                    return;

                ++_cur;

                if (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - _time_point).count() > delay || 
                    _cur == _n) {
                    write_progress();
                    _time_point = std::chrono::system_clock::now();
                }
            }

            const size_t& size() const {
                return _n;
            }

            const size_t& current() const {
                return _cur;
            }

            auto begin() {
                return iterator(this, _cur);
            }

            auto end() {
                return iterator(this, _n);
            }

            void set_description(const std::string& description) {
                _description = description;
                write();
            }

            void set_description(std::string&& description) {
                _description = std::move(description);
                write();
            }

            template<typename T>
            void set_description(const T& value) {
                set_description(boost::lexical_cast<std::string>(value));
            }

            void clear() const {
                if (_first)
                    return;

                printf("\033[2K\033[1F");
                fflush(stdout);
            }

        private:

            static constexpr size_t delay = 100;
            static constexpr auto long_dots = "................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................";
            static constexpr auto long_eqqs = "================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================";

#ifdef _WIN32
            static void _enable_vt100() {
                static bool enabled = false;
                if (!enabled) {
                    enabled = true;
                    system(" ");
                }
            }
#endif

            on_end _end;
            const size_t _n;
            std::string _description;
            bool _first = true, _enabled;
            size_t _cur = 0, _last = 0, _elapsed = 0, _k = 0;

            char* _buffer = new char[512];
            std::chrono::system_clock::time_point _time_point;


            static size_t get_window_size() {
#ifdef _WIN32
                CONSOLE_SCREEN_BUFFER_INFO csbi;
                GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
                return static_cast<size_t>(csbi.srWindow.Right) - csbi.srWindow.Left + 1;
#else
                struct winsize w;
                ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
                return w.ws_col;
#endif        
            }

            void write_progress() {
                if (!_enabled)
                    return;

                if (_first) {
                    printf("\n");
                    _first = false;
                    _time_point = std::chrono::system_clock::now();
                }

                const auto count = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - _time_point).count();
                const auto ips = count ? 1000 * (_cur - _last) / static_cast<double>(count): double(0);

                const auto est = count ? static_cast<size_t>((_n - _cur) / ips) : 0;
                const auto ht = est / 3600;
                const auto mt = (est - ht * 3600) / 60;
                const auto st = est - ht * 3600 - mt * 60;

                _elapsed += count;
                const auto es = _elapsed / 1000;
                const auto he = es / 3600;
                const auto me = (es - he * 3600) / 60;
                const auto se = es - he * 3600 - me * 60;

                _k = sprintf(_buffer, "| %zu/%zu [ %02zu:%02zu:%02zu<%02zu:%02zu:%02zu ] %.2f it/s ", _cur, _n, he, me, se, ht, mt, st, ips);

                write();

                _last = _cur;
            }

            void write() const {
                if (!_enabled)
                    return;
                
                const auto size = get_window_size();

                const auto progress = double(_cur) / _n;

                const auto has_description = static_cast<int>(_description.size() > 0);
                const auto used = _k + static_cast<int>(_description.size()) + has_description + 7;
                const auto middle = size > used ? size - used : 0;
                const auto eqqs = static_cast<int>(middle * progress);
                const auto dots = middle > eqqs ? middle - eqqs - 1 : 0;

                printf("\033[2K\033[1G%s%.*s %3d%%|%.*s%.*s%.*s%s",  _description.c_str(), has_description, ":",
                    static_cast<int>(progress * 100), eqqs, long_eqqs, static_cast<int>(_n != _cur), ">", static_cast<int>(dots), long_dots, _buffer);
                fflush(stdout);
            }

            class iterator {

            public:

                iterator(progress_bar* owner, const size_t& cur) : _owner(owner), _cur(cur) {}

                iterator& operator++() {
                    _owner->next();
                    ++_cur;
                    return *this;
                }

                bool operator==(const iterator& other) const {
                    return _owner == other._owner && _cur == other._cur;
                }

                bool operator!=(const iterator& other) const {
                    return !(*this == other);
                }

                const size_t& operator*() const {
                    return _cur;
                }

            private:

                size_t _cur;
                progress_bar* _owner;

            };

        };

    }// namespace utils

}// namespace acstc
