#pragma once

#include <utility>
#include <type_traits>
#include "fftw3.h"

namespace ample::utils {

    namespace _impl {

        template<bool Complex, bool Forward>
        struct plan;

        template<>
        struct plan<true, true> {

            template<typename T>
            static auto create(const int& size, T* data_in, T* data_out) {
                return fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(data_in),
                    reinterpret_cast<fftw_complex*>(data_out), FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
            }

        };

        template<>
        struct plan<true, false> {

            template<typename T>
            static auto create(const int& size, T* data_in, T* data_out) {
                return fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(data_out),
                    reinterpret_cast<fftw_complex*>(data_in), FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
            }

        };

        template<>
        struct plan<false, true> {

            template<typename T, typename V>
            static auto create(const int& size, T* data_in, V* data_out) {
                return fftw_plan_dft_r2c_1d(size, reinterpret_cast<double*>(data_in),
                    reinterpret_cast<fftw_complex*>(data_out), FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
            }

        };

        template<>
        struct plan<false, false> {

            template<typename T, typename V>
            static auto create(const int& size, T* data_in, V* data_out) {
                return fftw_plan_dft_c2r_1d(size, reinterpret_cast<fftw_complex*>(data_out),
                    reinterpret_cast<double*>(data_in), FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
            }

        };

        template<typename T, typename V, bool Complex>
        class basic_fft {

        private:

            using in_t = std::conditional_t<Complex, V, T>;
            using out_t = V;

        public:

            basic_fft() = delete;
            basic_fft(const basic_fft&) = delete;

            basic_fft(basic_fft&& other) noexcept {
                *this = std::move(other);
            }

            explicit basic_fft(const int& size) : _size(size), _dsize(size), _ssize(std::sqrt(_dsize)) {
                _forward_data = static_cast<in_t*>(fftw_malloc(sizeof(in_t) * size));
                _forward_data_end = _forward_data + size;

                _backward_data = static_cast<out_t*>(fftw_malloc(sizeof(out_t) * size));
                _backward_data_end = _backward_data + size;

                _plan_forward = plan<Complex, true>::create(size, _forward_data, _backward_data);
                _plan_backward = plan<Complex, false>::create(size, _forward_data, _backward_data);
            }

            ~basic_fft() {
                fftw_destroy_plan(_plan_forward);
                fftw_destroy_plan(_plan_backward);
                fftw_free(_forward_data);
                fftw_free(_backward_data);
            }

            basic_fft& operator=(basic_fft&& other) noexcept {
                std::swap(_size, other._size);
                std::swap(_dsize, other._dsize);
                std::swap(_forward_data, other._forward_data);
                std::swap(_backward_data, other._backward_data);
                std::swap(_forward_data_end, other._forward_data_end);
                std::swap(_backward_data_end, other._backward_data_end);
                std::swap(_plan_forward, other._plan_forward);
                std::swap(_plan_backward, other._plan_backward);
                return *this;
            }

            basic_fft& execute_forward() {
                fftw_execute(_plan_forward);
                return *this;
            }

            basic_fft& execute_backward() {
                fftw_execute(_plan_backward);
                return *this;
            }

            basic_fft& normalize_forward() {
                std::transform(_forward_data, _forward_data_end, _forward_data,
                   [this](const auto& value) { return value / _dsize; });
                return *this;
            }

            basic_fft& normalize_backward() {
                std::transform(_backward_data, _backward_data_end, _backward_data,
                   [this](const auto& value) { return value / _dsize; });
                return *this;
            }

            basic_fft& half_normalize_forward() {
                std::transform(_forward_data, _forward_data_end, _forward_data,
                   [this](const auto& value) { return value / _ssize; });
                return *this;
            }

            basic_fft& half_normalize_backward() {
                std::transform(_backward_data, _backward_data_end, _backward_data,
                   [this](const auto& value) { return value / _ssize; });
                return *this;
            }

            in_t* forward_data() {
                return _forward_data;
            }

            const in_t* forward_data() const {
                return _forward_data;
            }

            in_t& forward_data(const size_t& index) {
                return _forward_data[index];
            }

            const in_t& forward_data(const size_t& index) const {
                return _forward_data[index];
            }

            out_t* backward_data() {
                return _backward_data;
            }

            const out_t* backward_data() const {
                return _backward_data;
            }

            out_t& backward_data(const size_t& index) {
                return _backward_data[index];
            }

            const out_t& backward_data(const size_t& index) const {
                return _backward_data[index];
            }

            const in_t* forward_data_end() const {
                return _forward_data_end;
            }

            const out_t* backward_data_end() const {
                return _backward_data_end;
            }

            [[nodiscard]] const int& size() const {
                return _size;
            }

        private:

            int _size = 0;
            T _dsize, _ssize = 0;

            in_t* _forward_data = nullptr, *_forward_data_end = nullptr;
            out_t* _backward_data = nullptr, *_backward_data_end = nullptr;
            fftw_plan _plan_forward = fftw_plan(), _plan_backward = fftw_plan();

        };


    }// namespace _impl


    template<typename T, typename V>
    using real_fft = _impl::basic_fft<T, V, false>;

    template<typename T, typename V>
    using complex_fft = _impl::basic_fft<T, V, true>;

}// namespace ample::utils
