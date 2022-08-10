#pragma once
#include <tuple>
#include <variant>
#include <functional>
#include "utils/assert.hpp"

namespace ample::threads {

    template<typename...>
    class pool;

    template<typename... Parameters>
    class task {

    private:

        using owner_t = pool<Parameters...>;

    public:

        [[nodiscard]] size_t id() const {
            return _id;
        }

        void push(typename owner_t::template parameter_value<Parameters>... parameters) const {
            _owner.push(_id, std::move(parameters)...);
        }

        void push_single(typename owner_t::template parameter_value<Parameters>... parameters) const {
            _owner.push_single(_id, std::move(parameters)...);
        }

        void remove() const {
            _owner.remove(_id);
        }

        void wait() const {
            _owner.wait(_id);
        }

        [[nodiscard]] owner_t& owner() {
            return _owner;
        }

        [[nodiscard]] const owner_t& owner() const {
            return _owner;
        }

    private:

        template<typename F>
        task(owner_t& owner, const size_t& id, F function) : _id(id), _owner(owner), _function(std::forward<F>(function)) {}

        friend owner_t;

        size_t _id;
        owner_t& _owner;
        std::variant<std::function<void(Parameters&...)>, std::function<void(Parameters&..., const task&)>> _function;

        void _execute(Parameters&... parameters) {
            switch (_function.index()) {
                case 0:
                    std::get<std::function<void(Parameters&...)>>(_function)(parameters...);
                    break;
                case 1:
                    std::get<std::function<void(Parameters&..., const task&)>>(_function)(parameters..., *this);
                    break;
                default:
                    utils::dynamic_assert(false, "Task function variant does not hold a value");

            }
        }

    };

}// namespace ample::threads
