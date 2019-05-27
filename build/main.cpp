#include <cmath>
#include <chrono>
#include <string>
#include <cstddef>
#include <iostream>
#include <algorithm>
#include "config.hpp"
#include "solver.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/callback.hpp"
#include "initial_conditions.hpp"
#include "boost/program_options.hpp"

namespace types = acstc::types;
namespace po = boost::program_options;

template<typename KS, typename PS>
auto get_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s) {
    KS ws(k0.size());
    PS as(phi_s.size());
    std::transform(phi_s.begin(), phi_s.end(), as.begin(), [](const auto& phi) { return phi / (2 * std::sqrt(M_PI)); });
    std::transform(k0.begin(), k0.end(), ws.begin(), [](const auto& k0) { return 1 / std::pow(k0, 2); } );
    return acstc::greene_source<types::real_t>(config.y0(), config.y1(), config.ny(), config.y_s(), as, ws);
}

template<typename F>
auto add_verbosity(const size_t& verbosity, const size_t& report, F& function) {
    return [&](const auto& k_j, const auto& phi_j, auto&& writer) mutable {
        if (verbosity) {
            const auto start = std::chrono::system_clock::now();
            if (report)
                function(k_j, phi_j, acstc::utils::callbacks(writer, acstc::utils::progress_callback(report)));
            else
                function(k_j, phi_j, writer);
            const auto end = std::chrono::system_clock::now();
            std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << "ms" << std::endl;
            return;
        }
        function(k_j, phi_j, writer);
    };
}

template<typename F>
auto add_writer(const acstc::config<types::real_t>& config, const std::string& filename,
        const bool binary, const size_t& step, F& function) {
    return [&, binary](const auto& k_j, const auto& phi_j) mutable {
        const auto xs = acstc::utils::mesh_1d(config.x0(), config.x1(), (config.nx() - 1) / step + 1);
        const auto ys = acstc::utils::mesh_1d(config.y0(), config.y1(), config.ny());
        if (binary) {
            const auto nx = static_cast<const uint32_t>(xs.size());
            const auto ny = static_cast<const uint32_t>(ys.size());
            acstc::utils::binary_writer<double> writer(filename);
            writer.stream().write(reinterpret_cast<const char*>(&nx), sizeof(nx));
            writer.stream().write(reinterpret_cast<const char*>(&ny), sizeof(ny));
            writer.write(xs);
            writer.write(ys);
            function(k_j, phi_j, acstc::utils::ekc_callback(step, [&](const auto data) {
                writer.write(reinterpret_cast<const double*>(data.data()), ny * 2);
            }));
            return;
        }
        acstc::utils::text_writer<types::real_t> writer(filename);
        writer.write(xs);
        writer.write(ys);
        function(k_j, phi_j, acstc::utils::ekc_callback(step, [&](const auto data) {
            writer.write(reinterpret_cast<const double*>(data.data()), ys.size() * 2);
        }));
    };
}

template<typename T, typename F>
auto add_modes(const acstc::config<types::real_t>& config, F& function) {
    return [&]() mutable {
        if (config.const_modes()) {
            const auto [k_j, phi_j] = config.create_const_modes<T>();
            function(k_j, phi_j);
            return;
        }
        const auto [k_j, phi_j] = config.create_modes<T>();
        function(k_j, phi_j);
    };
}

int main(int argc, char* argv[]) {
    po::positional_options_description positional;
    positional.add("domain_config_filename", 1);

    po::options_description generic("Generic options");
    size_t verbosity, report;
    generic.add_options()
            ("help,h", "Print this message")
            ("verbosity,v", po::value(&verbosity)->default_value(0), "Verbosity level (0-1)")
            ("report,r", po::value(&report)->default_value(0)->value_name("k"),
                    "If verbosity level > 0 report every k computed rows (0 = don't report)");

    po::options_description output("Output options");
    size_t step;
    std::string output_filename;
    output.add_options()
            ("output,o", po::value(&output_filename)->default_value("output.txt"), "Output filename")
            ("step,s", po::value(&step)->default_value(100)->value_name("k"), "Output every k-th computed row")
            ("binary", "Use binary output");

    po::options_description options;
    std::string config_filename;
    options.add_options()("domain_config_filename", po::value(&config_filename)->default_value("config.json"));
    options.add(generic).add(output);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).positional(positional).options(options).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        po::options_description desc;
        desc.add(generic).add(output);
        std::cout << "Usage: domain_config_filename (=config.json) [options]\n" << desc << std::endl;
        return 0;
    }

    acstc::config config(config_filename);
    acstc::solver solver(config);

    const auto k0 = config.k0();
    const auto phi_s = config.phi_s();

    const auto init = get_initial_conditions(config, k0, phi_s);

    auto with_solver = [&](const auto& k_j, const auto& phi_j, auto&& callback) mutable {
        solver.solve(init, k0, k_j, phi_j, callback);
    };

    auto with_verbosity = add_verbosity(verbosity, report, with_solver);
    auto with_writer = add_writer(config, output_filename, vm.count("binary") > 0, step, with_verbosity);

    if (config.complex_modes())
        add_modes<types::complex_t>(config, with_writer)();
    else
        add_modes<types::real_t>(config, with_writer)();
}