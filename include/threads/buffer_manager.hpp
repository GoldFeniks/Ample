#pragma once
#include <mutex>
#include <shared_mutex>
#include <condition_variable>
#include "utils/types.hpp"

namespace ample::threads {

    template<typename T>
    class buffer_manager {

    private:

        class buffer_wrapper {

        public:

            buffer_wrapper(const buffer_wrapper&) = delete;

            virtual ~buffer_wrapper() {
                if (_owner != nullptr)
                    unlock();
            }

            void unlock() {
                if (_lock)
                    _lock.unlock();
            }

            void lock() {
                if (!_lock)
                    _lock.lock();
            }

            T& data() {
                return _owner->_data[_buffer_index];
            }

            const T& data() const {
                return _owner->_data[_buffer_index];
            }

            [[nodiscard]] bool ready() const {
                return std::all_of(_owner->_done[_buffer_index].begin(), _owner->_done[_buffer_index].end(), [](const auto& v) { return v; });
            }

            virtual bool complete() = 0;

        protected:

            friend class buffer_manager;

            buffer_manager* _owner = nullptr;
            std::unique_lock<std::mutex> _lock;
            size_t _buffer_index;

            buffer_wrapper(buffer_manager* owner, const size_t& buffer_index) : _owner(owner), _lock(owner->_mutex[buffer_index]), _buffer_index(buffer_index) {}

        };

    public:

        class thread_buffer_wrapper : public buffer_wrapper {

        public:

            using buffer_wrapper::ready;

            bool complete() override {
                _owner->_done[_buffer_index][_thread_index] = true;
                _owner->_index[_thread_index] = (_owner->_index[_thread_index] + 1) % _owner->size();

                return ready();
            }

            [[nodiscard]] bool available() const {
                return !_owner->_done[_buffer_index][_thread_index];
            }

            void await() {
                _owner->_cv[_buffer_index].wait(_lock, [this]{ return available(); });
            }

        private:

            friend class buffer_manager;

            thread_buffer_wrapper(buffer_manager* owner, const size_t& buffer_index, const size_t& thread_index) :
                buffer_wrapper(owner, buffer_index), _thread_index(thread_index) {}

            size_t _thread_index;
            using buffer_wrapper::_lock;
            using buffer_wrapper::_owner;
            using buffer_wrapper::_buffer_index;

        };

        class result_buffer_wrapper : public buffer_wrapper {

        public:

            using buffer_wrapper::ready;

            ~result_buffer_wrapper() {
                if (_owner != nullptr)
                    _owner->_cv[_buffer_index].notify_all();
            }

            bool complete() override {
                _owner->_done[_buffer_index].assign(_owner->_done[_buffer_index].size(), false);

                std::lock_guard<std::mutex> lock(_owner->_index_mutex);
                _owner->_current_index = (_owner->_current_index + 1) % _owner->size();

                return false;
            }

        private:

            friend class buffer_manager;

            result_buffer_wrapper(buffer_manager* owner, const size_t& buffer_index) : buffer_wrapper(owner, buffer_index) {}

            using buffer_wrapper::_owner;
            using buffer_wrapper::_buffer_index;

        };

        buffer_manager(const size_t& n_workers, const size_t& buffer_size, const T& initial_value) :
            _data(buffer_size, initial_value),
            _done(buffer_size, types::vector1d_t<bool>(n_workers, false)),
            _index(n_workers, 0),
            _mutex(buffer_size),
            _cv(buffer_size) {}

        thread_buffer_wrapper get(const size_t& i) {
            return thread_buffer_wrapper(this, _index[i], i);
        }

        result_buffer_wrapper current() {
            std::lock_guard<std::mutex> index_lock(_index_mutex);
            return result_buffer_wrapper(this, _current_index);
        }

        [[nodiscard]] size_t size() const {
            return _data.size();
        }

    private:

        std::mutex _index_mutex;
        size_t _current_index = 0;
        types::vector1d_t<T> _data;
        types::vector2d_t<bool> _done;
        types::vector1d_t<size_t> _index;
        types::vector1d_t<std::mutex> _mutex;
        types::vector1d_t<std::condition_variable> _cv;

    };

}// namespace ample::threads
