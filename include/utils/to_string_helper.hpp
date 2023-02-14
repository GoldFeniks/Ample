#pragma once

namespace ample::utils {

    class to_string_helper {

    public:

        to_string_helper() {
            _stream.setf(std::ios_base::boolalpha);
            _stream.unsetf(std::ios_base::showpoint);
        }

        template<typename V>
        std::string to_string(const V& value) {
            _stream << value;
            const auto result = _stream.str();
            _stream.str(std::string());
            return result;
        }

    private:

        std::stringstream _stream{};

    };

}// namespace ample::utils
