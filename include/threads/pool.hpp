#pragma once
#include <queue>
#include <mutex>
#include <thread>
#include <optional>
#include <type_traits>
#include <unordered_map>
#include <condition_variable>
#include "task.hpp"
#include "utils/types.hpp"
#include "utils/assert.hpp"

namespace ample::threads {

    static constexpr struct resource_t {} resource{};

    template<typename... Parameters>
    class pool {

    private:

        using task_t = task<Parameters...>;
        using parameters_t = std::tuple<Parameters...>;

        friend task_t;

        template<typename T>
        struct parameter_value {

            std::variant<T, resource_t> value;

            parameter_value() = default;
            parameter_value(T value) : value(std::move(value)) {}
            parameter_value(resource_t value) : value(value) {}

            T& get(std::optional<types::vector1d_t<T>>& res, const size_t& index) {
                switch (this->value.index()) {
                    case 0:
                        return std::get<T>(value);
                    case 1:
                        return res.value()[index];
                    default:
                        throw std::runtime_error("No parameter value");
                }
            }

        };

    public:

        static constexpr resource_t resource{};

        explicit pool(const size_t& workers) {
            utils::dynamic_assert(workers > 0, "Number of workers must be greater than 0");

            for (size_t i = 0; i < workers; ++i)
                _workers.emplace_back(&pool::_worker, this, i);
        }

        template<typename F>
        task_t& add(F function) {
            std::lock_guard lock(_mutex);
            auto [it, added] = _tasks.try_emplace(_id, *this, _id, std::forward<F>(function));
            ++_id;
            return it->second.task;
        }

        template<typename F, typename... P>
        task_t& execute(F function, P... parameters) {
            auto& task = add(std::forward<F>(function));
            push(task._id, std::forward<P>(parameters)...);
            return task;
        }

        [[nodiscard]] size_t workers() const {
            return _workers.size();
        }

        void join() {
            {
                std::lock_guard lock(_mutex);
                _stop = true;
            }

            for (auto& it : _workers)
                it.join();
        }

        void push(const size_t& id, parameter_value<Parameters>... parameters) {
            _push(id, false, std::move(parameters)...);
        }

        void push(const task_t& task, parameter_value<Parameters>... parameters) {
            push(task._id, false, std::move(parameters)...);
        }

        void push_single(const size_t& id, parameter_value<Parameters>... parameters) {
            _push(id, true, std::move(parameters)...);
        }

        void push_single(const task_t& task, parameter_value<Parameters>... parameters) {
            push(task._id, true, std::move(parameters)...);
        }

        [[nodiscard]] task_t& get_task(const size_t& id) {
            return _find(id).task;
        }

        [[nodiscard]] const task_t& get_task(const size_t& id) const {
            return _find(id).task;
        }

        [[nodiscard]] bool empty() const {
            std::lock_guard lock(_mutex);
            return _queue.empty();
        }

        void remove(const size_t& id) {
            std::lock_guard lock(_mutex);
            auto& holder = _find(id);
            holder.remove = true;
            _check_remove(holder);
        }

        void remove(const task_t& task) {
            remove(task._id);
        }

        [[nodiscard]] bool contains(const size_t& id) {
            return _tasks.contains(id);
        }

        void wait(const size_t& id) {
            std::unique_lock lock(_mutex);
            auto& holder = _find(id);

            _completion_cv.wait(lock, [&, this]{ return !contains(id) || holder.usages == 0; });
        }

        template<size_t I>
        void add_resource(std::tuple_element_t<I, parameters_t> value) {
            std::get<I>(_resources) = std::vector<std::tuple_element_t<I, parameters_t>>(workers(), value);
        }

        void pause() {
            std::lock_guard lock(_mutex);
            _pause = true;
        }

        void resume() {
            std::lock_guard lock(_mutex);
            _pause = false;
            _cv.notify_all();
        }

    private:

        struct task_holder {

            task_t task;
            size_t usages{};
            bool remove = false;

            template<typename F>
            task_holder(pool& owner, const size_t& id, F function) : task(task_t(owner, id, std::forward<F>(function))) {}

        };

        size_t _id{};
        types::vector1d_t<std::thread> _workers;
        std::unordered_map<size_t, task_holder> _tasks;
        std::tuple<std::optional<types::vector1d_t<Parameters>>...> _resources;
        std::queue<std::pair<size_t, std::tuple<parameter_value<Parameters>...>>> _queue;

        bool _stop = false;
        bool _pause = false;
        std::mutex _mutex;
        std::condition_variable _cv, _completion_cv;

        void _worker(size_t i) {
            task_holder* holder = nullptr;
            std::pair<size_t, std::tuple<parameter_value<Parameters>...>> info;

            while (true) {
                {
                    std::unique_lock lock(_mutex);
                    _cv.wait(lock, [this] { return !_queue.empty() && !_pause || _stop; });

                    if (_queue.empty() && _stop)
                        break;

                    info = std::move(_queue.front());
                    _queue.pop();

                    holder = &_find(info.first);
                }

                _apply(holder->task, info.second, i, std::make_index_sequence<sizeof...(Parameters)>{});

                {
                    std::unique_lock lock(_mutex);

                    if (!--holder->usages)
                        _completion_cv.notify_all();

                    _check_remove(*holder);
                }

                holder = nullptr;
            }
        }

        task_holder& _find(const size_t& id) {
            const auto it = _tasks.find(id);
            utils::dynamic_assert(it != _tasks.end(), "Task with id ", id, " does not exist");
            return it->second;
        }

        const task_holder& _find(const size_t& id) const {
            const auto it = _tasks.find(id);
            utils::dynamic_assert(it != _tasks.end(), "Task with id ", id, " does not exist");
            return it->second;
        }

        void _check_remove(task_holder& holder) {
            if (holder.usages == 0 && holder.remove)
                _tasks.erase(holder.task._id);
        }

        template<typename... P>
        void _push(const size_t& id, const bool single, P... parameters) {
            std::lock_guard lock(_mutex);

            auto& holder = _find(id);

            if (single && holder.usages)
                return;

            ++holder.usages;

            _queue.emplace(id, std::make_tuple(parameter_value<Parameters>(std::forward<P>(parameters))...));

            if (!_pause)
                _cv.notify_one();
        }

        template<size_t... I>
        void _apply(task_t& task, std::tuple<parameter_value<Parameters>...>& parameters, const size_t& i, std::index_sequence<I...>) {
            task._execute(std::get<I>(parameters).get(std::get<I>(_resources), i)...);
        }

    };

}// namespace ample::threads
