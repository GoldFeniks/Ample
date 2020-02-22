#pragma once
#include <string>
#include <cstddef>
#include <iostream>
#include <algorithm>

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

            progress_bar(const size_t& n, std::string title) : 
                _n(n),
                _n_string(std::string("/") + std::to_string(n)),
                _title(std::move(title))
            {
                const auto size = get_window_size();
                _clear_string = std::string(size, ' ');

                _iterations_width = static_cast<int>(_n_string.size() - 1);

                const auto available = size - _n_string.size() * 2 - _title.size() - 9;
                _progress_string = std::string(available, '.');
            }

            void operator()() {
                next();
            }

            void next() {
                if (_cur >= _n)
                    return;

                ++_cur;

                if (_cur == _n) {
                    for (size_t i = _last_ind; i < _progress_string.size(); ++i)
                        _progress_string[i] = '=';
                    write_bar(100);
                    printf("\n");
                    return;
                }

                write_progress();
            }

            const size_t& size() const {
                return _n;
            }

            const size_t& current() const {
                return _cur;
            }

        private:

            size_t _cur = 0, _last_ind = 0;
            int _iterations_width;
            const size_t _n;
            std::string _n_string, _clear_string, _progress_string, _title;

            static size_t get_window_size() {
#ifdef _WIN32
                CONSOLE_SCREEN_BUFFER_INFO csbi;
                GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
                return csbi.srWindow.Right - csbi.srWindow.Left + 1;
#else
                struct winsize w;
                ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
                return w.ws_col;
#endif        
            }

            void write_progress() {
                const auto progress = double(_cur) / _n;
                const auto ind = size_t(progress * _progress_string.size());

                if (ind > _last_ind) {
                    for (size_t i = _last_ind; i < ind; ++i)
                        _progress_string[i] = '=';
                    _last_ind = ind;
                }
                _progress_string[ind] = '>';
                write_bar(progress * 100);
            }

            void write_bar(const size_t& progress) const {
                printf("\r%s\r%s: %3d%%|%s| %*zu%s ", _clear_string.c_str(), _title.c_str(), progress, _progress_string.c_str(), 
                    _iterations_width, _cur, _n_string.c_str());
                fflush(stdout);
            }

        };

    }// namespace utils


}// namespace acstc