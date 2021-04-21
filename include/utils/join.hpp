#pragma once
#include <string>
#include <sstream>

namespace ample::utils {

    template<typename It>
    std::string join_it(It begin, const It& end, const std::string& sep = " ") {
        std::ostringstream result;

        if (begin != end)
            result << *begin++;

        while (begin != end)
            result << sep << *begin++;

        return result.str();
    }

    template<typename First, typename... Tail>
    std::string join(const First& first, const Tail&... tail) {
        std::ostringstream result;
        ((result << first), ..., (result << tail));
        return result.str();
    }

}// namespace ample::utils