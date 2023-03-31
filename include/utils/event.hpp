#pragma once
#include <memory>
#include <functional>
#include <unordered_set>

namespace ample::utils {

    template<typename... Args>
    class event {

    public:

        using func_t = std::function<void(Args&...)>;

        std::weak_ptr<func_t> operator+=(func_t func) {
            const auto [pointer, _] = _callbacks.emplace(std::make_shared<func_t>(std::move(func)));
            return *pointer;;
        }

        void operator-=(const std::weak_ptr<func_t>& func) {
            _callbacks.erase(func.lock());
        }

        void operator()(Args... args) const {
            for (auto& it : _callbacks)
                it->operator()(args...);
        }

    private:

        mutable std::unordered_set<std::shared_ptr<func_t>> _callbacks;

    };

}// namespace ample::utils
