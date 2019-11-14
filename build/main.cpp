#include <cmath>
#include <chrono>
#include <string>
#include <cstddef>
#include <iostream>
#include <algorithm>
#include "config.hpp"
#include "solver.hpp"
#include "io/writer.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/callback.hpp"
#include "utils/verbosity.hpp"
#include "initial_conditions.hpp"
#include "utils/interpolation.hpp"
#include "boost/program_options.hpp"

namespace types = acstc::types;
namespace po = boost::program_options;

using acstc::utils::verboseln;

template<typename KS, typename PS>
auto get_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s) {
    KS ws(k0.size());
    PS as(phi_s.size());
    std::transform(phi_s.begin(), phi_s.end(), as.begin(), [](const auto& phi) { return phi / (2 * std::sqrt(M_PI)); });
    std::transform(k0.begin(), k0.end(), ws.begin(), [](const auto& k0) { return 1 / std::pow(k0, 2); } );
    return acstc::greene_source<types::real_t>(config.y0(), config.y1(), config.ny(), config.y_s(), as, ws);
}

template<typename F>
auto add_verbosity(const size_t& report, F& function) {
    return [&](auto&& writer) mutable {
        const auto start = std::chrono::system_clock::now();
        if (report)
            function(acstc::utils::callbacks(writer, acstc::utils::progress_callback(report)));
        else
            function(writer);
        const auto end = std::chrono::system_clock::now();
        verboseln(1, "Elapsed time: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");
        return;
    };
}

template<typename F>
auto add_writer(const acstc::config<types::real_t>& config, const std::string& filename,
        const bool binary, const size_t& step, F& function) {
    return [&, binary]() mutable {
        const auto xs = acstc::utils::mesh_1d(config.x0(), config.x1(), (config.nx() - 1) / step + 1);
        const auto ys = acstc::utils::mesh_1d(config.y0(), config.y1(), config.ny());
        if (binary) {
            const auto nx = static_cast<const uint32_t>(xs.size());
            const auto ny = static_cast<const uint32_t>(ys.size());
            acstc::utils::binary_writer<types::real_t> writer(filename);
            writer.stream().write(reinterpret_cast<const char*>(&nx), sizeof(nx));
            writer.stream().write(reinterpret_cast<const char*>(&ny), sizeof(ny));
            writer.write(xs);
            writer.write(ys);
            function(acstc::utils::ekc_callback(step, [&](const auto data) {
                writer.write(reinterpret_cast<const types::real_t*>(data.data()), ny * 2);
            }));
            return;
        }
        acstc::utils::text_writer<types::real_t> writer(filename);
        writer.stream() << xs.size() << ' ' << ys.size() << '\n';
        writer.write(xs);
        writer.write(ys);
        function(acstc::utils::ekc_callback(step, [&](const auto data) {
            writer.write(reinterpret_cast<const types::real_t*>(data.data()), ys.size() * 2);
        }));
    };
}

template<typename T, typename F1, typename F2>
auto add_modes(const acstc::config<types::real_t>& config, const size_t mn, F1& function, F2& function_const) {
    return [&, mn](auto&& callback) mutable {
        if (config.const_modes()) {
            const auto start = std::chrono::system_clock::now();
            auto [k_j, phi_j] = config.create_const_modes<T>();
            const auto end = std::chrono::system_clock::now();
            verboseln(1, "Modes computing time: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");

            if (k_j.size() > mn) {
                k_j.erase_last(k_j.size() - mn);
                phi_j.erase_last(phi_j.size() - mn);
            }
            function_const(k_j, phi_j, callback);
            return;
        }
        const auto start = std::chrono::system_clock::now();
        auto [k_j, phi_j] = config.create_modes<T>();
        const auto end = std::chrono::system_clock::now();
        verboseln(1, "Modes computing time: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");
        if (k_j.size() > mn) {
            k_j.erase_last(k_j.size() - mn);
            phi_j.erase_last(phi_j.size() - mn);
        }
        function(k_j, phi_j, callback);
    };
}

template<typename V, typename W>
void write_modes(const acstc::utils::linear_interpolated_data_1d<types::real_t, V>& modes, W& writer) {
    for (size_t i = 0; i < modes.size(); ++i)
        writer(modes[i].data());
}

template<typename V, typename W>
void write_modes(const acstc::utils::linear_interpolated_data_2d<types::real_t, V>& modes, W& writer) {
    for (size_t i = 0; i < modes.size(); ++i)
        for (const auto& it : modes[i].data())
            writer(it);
}

template<typename KV, typename PV, typename W, typename WR>
void write_modes(const KV& k_j, const PV& phi_j, W& writer, const WR& wrapper) {
    auto wrapped = wrapper(writer);
    write_modes(k_j, wrapped);
    write_modes(phi_j, writer);
}

template<typename T, typename W>
void save_modes(const acstc::config<types::real_t>& config, const std::string& filename, const bool binary, const W& wrapper) {
    if (config.const_modes()) {
        const auto [k_j, phi_j] = config.create_const_modes<T>();
        const auto& ys = k_j.template get<0>();
        if (binary) {
            const auto ny = static_cast<const uint32_t>(ys.size());
            const auto nm = static_cast<const uint32_t>(k_j.size());
            acstc::utils::binary_writer<types::real_t> writer(filename);
            writer.stream().write(reinterpret_cast<const char*>(&ny), sizeof(ny));
            writer.stream().write(reinterpret_cast<const char*>(&nm), sizeof(nm));
            writer.write(ys);
            write_modes(k_j, phi_j, writer, wrapper);
            return;
        }
        acstc::utils::text_writer<types::real_t> writer(filename);
        writer.stream() << ys.size() << ' ' << k_j.size() << '\n';
        writer.write(ys);
        write_modes(k_j, phi_j, writer, wrapper);
        return;
    }
    const auto [k_j, phi_j] = config.create_modes<T>();
    const auto& xs = k_j.template get<0>();
    const auto& ys = k_j.template get<1>();
    if (binary) {
            const auto nx = static_cast<const uint32_t>(xs.size());
            const auto ny = static_cast<const uint32_t>(ys.size());
            const auto nm = static_cast<const uint32_t>(k_j.size());
            acstc::utils::binary_writer<types::real_t> writer(filename);
            writer.stream().write(reinterpret_cast<const char*>(&nx), sizeof(nx));
            writer.stream().write(reinterpret_cast<const char*>(&ny), sizeof(ny));
            writer.stream().write(reinterpret_cast<const char*>(&nm), sizeof(nm));
            writer.write(xs);
            writer.write(ys);
            write_modes(k_j, phi_j, writer, wrapper);
            return;
        }
    acstc::utils::text_writer<types::real_t> writer(filename);
    writer.stream() << xs.size() << ' ' << ys.size() << ' ' << k_j.size() << '\n';
    writer.write(xs);
    writer.write(ys);
    write_modes(k_j, phi_j, writer, wrapper);
}

int main(int argc, char* argv[]) {
    po::positional_options_description positional;
    positional.add("job_type", 1);

    po::options_description generic("Generic options");
    size_t report;
    std::string config_filename;
    generic.add_options()
            ("help,h", "Print this message")
            ("verbosity,v", po::value(&acstc::utils::verbosity::instance().level)->default_value(0), "Verbosity level")
            ("report,r", po::value(&report)->default_value(0)->value_name("k"),
                    "If verbosity level > 0 report every k computed rows (0 = don't report)")
            ("config,c", po::value(&config_filename)->default_value("config.json"), "Config filename");

    po::options_description output("Output options");
    size_t step;
    std::string output_filename;
    output.add_options()
            ("output,o", po::value(&output_filename)->default_value("output.txt"), "Output filename")
            ("step,s", po::value(&step)->default_value(100)->value_name("k"), "Output every k-th computed row")
            ("binary", "Use binary output");

    po::options_description calculation("Calculation options");
    size_t num_workers, buff_size;
    calculation.add_options()
            ("workers,w", po::value(&num_workers)->default_value(1), "Number of workers for calculation")
            ("buff,b", po::value(&buff_size)->default_value(100), "Buff size to be used during multithreaded calculation");

    po::options_description options;
    std::string job_type;
    options.add_options()
            ("job_type", po::value(&job_type)->default_value("solution"));
    options.add(generic).add(output).add(calculation);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).positional(positional).options(options).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        po::options_description desc;
        desc.add(generic).add(output).add(calculation);
        std::cout << "Usage: [solution|modes] (=solution) [options]\n" << desc << std::endl;
        return 0;
    }

    acstc::config config(config_filename);

    if (job_type == "solution") {
        acstc::solver solver(config);

        const auto [k0, phi_s] = config.create_source_modes();

        const auto init = get_initial_conditions(config, k0, phi_s);

        auto execute_function = [&](auto&& with_solver) {
            auto with_verbosity = add_verbosity(report, with_solver);
            auto with_writer = add_writer(config, output_filename, vm.count("binary") > 0, step, with_verbosity);
            with_writer();
        };

        auto with_solver = [&](const auto& k_j, const auto& phi_j, auto&& callback) mutable {
            solver.solve(init, k0, k_j, phi_j, callback, config.border_width(), config.past_n(), num_workers, buff_size);
        };

        auto with_solver_const = [&](const auto& k_j, const auto& phi_j, auto&& callback) mutable {
            solver.solve(init, k0, k_j, phi_j, callback, config.past_n(), num_workers, buff_size);
        };

        if (config.complex_modes())
            execute_function(add_modes<types::complex_t>(config, k0.size(), with_solver, with_solver_const));
        else
            execute_function(add_modes<types::real_t>(config, k0.size(), with_solver, with_solver_const));
        return 0;
    }
    if (job_type == "modes") {
        if (config.complex_modes())
            save_modes<types::complex_t>(config, output_filename, vm.count("binary") > 0,
                [](auto& writer) {
                    return [&](const auto& data) {
                        writer.write(reinterpret_cast<const types::real_t*>(data.data()), data.size() * 2);
                    };
            });
        else
            save_modes<types::real_t>(config, output_filename, vm.count("binary") > 0,
                    [](auto& writer) { return [&](const auto& data) { writer(data); }; });
        return 0;
    }
    throw std::logic_error(std::string("Unknown job type: ") + job_type);
}
