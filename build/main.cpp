#include <cmath>
#include <chrono>
#include <string>
#include <cstddef>
#include <iostream>
#include <algorithm>
#include "rays.hpp"
#include "config.hpp"
#include "solver.hpp"
#include "io/writer.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "utils/callback.hpp"
#include "utils/verbosity.hpp"
#include "initial_conditions.hpp"
#include "utils/progress_bar.hpp"
#include "utils/interpolation.hpp"
#include "boost/program_options.hpp"

namespace types = acstc::types;
namespace po = boost::program_options;

using acstc::utils::verboseln;

template<typename KS, typename PS>
auto get_ray_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s,
    const acstc::utils::linear_interpolated_data_1d<types::real_t>& k_j) {
    return acstc::ray_source(config.x0(), config.y0(), config.y1(), config.ny(), 0., config.y_s(), config.l1(), config.nl(),
            config.a0(), config.a1(), config.na(), k0, phi_s, k_j);
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s) {
    auto [k_j, phi_j] = config.create_const_modes<types::real_t>();

    if (k_j.size() > k0.size())
        k_j.erase_last(k_j.size() - k0.size());

    return get_ray_initial_conditions(config, k0, phi_s, phi_j);
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s,
    const acstc::utils::linear_interpolated_data_1d<types::real_t, types::complex_t>& k_j) {
    const auto& ys = k_j.template get<0>();
    types::vector2d_t<types::real_t> new_k_j(k_j.size(), types::vector1d_t<types::real_t>(ys.size()));

    for (size_t j = 0; j < k_j.size(); ++j)
        std::transform(k_j[j].data().begin(), k_j[j].data().end(), new_k_j[j].begin(), [](const auto& v) { return v.real(); });

    return get_ray_initial_conditions(config, k0, phi_s, acstc::utils::linear_interpolated_data_1d<types::real_t>(ys, new_k_j));
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s,
    const acstc::utils::linear_interpolated_data_2d<types::real_t>& k_j) {
    types::vector2d_t<types::real_t> new_k_j;
    new_k_j.reserve(k_j.size());

    for (size_t j = 0; j < k_j.size(); ++j)
        new_k_j.emplace_back(k_j[j][0].begin(), k_j[j][0].end());

    return get_ray_initial_conditions(config, k0, phi_s, 
        acstc::utils::linear_interpolated_data_1d<types::real_t>(k_j.template get<1>(), new_k_j));
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s,
    const acstc::utils::linear_interpolated_data_2d<types::real_t, types::complex_t>& k_j) {
    const auto& ys = k_j.template get<1>();
    types::vector2d_t<types::real_t> new_k_j(k_j.size(), types::vector1d_t<types::real_t>(ys.size()));

    for (size_t j = 0; j < k_j.size(); ++j)
        std::transform(k_j[j][0].begin(), k_j[j][0].end(), new_k_j[j].begin(), [](const auto& v) { return v.real(); });

    return get_ray_initial_conditions(config, k0, phi_s, acstc::utils::linear_interpolated_data_1d<types::real_t>(ys, new_k_j));
}

template<typename KS, typename PS>
auto get_simple_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s) {
    const auto init = config.init();

    KS ws(k0.size());
    PS as(phi_s.size());
    std::transform(phi_s.begin(), phi_s.end(), as.begin(), [](const auto& phi) { return phi / (2 * std::sqrt(M_PI)); });
    std::transform(k0.begin(), k0.end(), ws.begin(), [](const auto& k0) { return 1 / std::pow(k0, 2); } );

    if (init == "green")
        return acstc::greene_source<types::complex_t>(config.y0(), config.y1(), config.ny(), config.y_s(), as, ws);

    if (init == "gauss")
        return acstc::gaussian_source<types::complex_t>(config.y0(), config.y1(), config.ny(), config.y_s(), as, ws);

    throw new std::logic_error(std::string("Unknown initial conditions type: ") + init);
}

template<typename KS, typename PS, typename KJ>
auto get_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s, const KJ& k_j) {
    const auto init = config.init();

    if (init == "ray")
        return get_ray_initial_conditions(config, k0, phi_s, k_j);

    return get_simple_initial_conditions(config, k0, phi_s);
}

template<typename KS, typename PS>
auto get_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s) {
    const auto init = config.init();

    if (init == "ray")
        return get_ray_initial_conditions(config, k0, phi_s);

    return get_simple_initial_conditions(config, k0, phi_s);
}

template<typename F>
auto add_verbosity(const size_t& report, const size_t& n, F& function) {
    return [&](auto&& writer) mutable {
        const auto start = std::chrono::system_clock::now();
        if (acstc::utils::verbosity::instance().level >= 2)
            function(acstc::utils::callbacks(writer, acstc::utils::progress_bar_callback(n, "Solution")));
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
            auto [k_j, phi_j] = config.create_const_modes<T>(acstc::utils::verbosity::instance().level >= 2);
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
        auto [k_j, phi_j] = config.create_modes<T>(acstc::utils::verbosity::instance().level >= 2);
        const auto end = std::chrono::system_clock::now();
        verboseln(1, "Modes computing time: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");
        if (k_j.size() > mn) {
            k_j.erase_last(k_j.size() - mn);
            phi_j.erase_last(phi_j.size() - mn);
        }
        function(k_j, phi_j, callback);
    };
}

template<typename KS, typename PS, typename F>
auto add_initial_conditions(const acstc::config<types::real_t>& config, const KS& k0, const PS& phi_s, F& function) {
    return [&](const auto& k_j, const auto& phi_j, auto&& callback) mutable {
        function(get_initial_conditions(config, k0, phi_s, k_j), k_j, phi_j, callback);
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
        const auto [k_j, phi_j] = config.create_const_modes<T>(acstc::utils::verbosity::instance().level >= 2);
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
    const auto [k_j, phi_j] = config.create_modes<T>(acstc::utils::verbosity::instance().level >= 2);
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

template<typename V, typename W>
void write_strided(const V& v, const size_t& k, W& writer) {
    auto [begin, end] = acstc::utils::stride(v.begin(), v.end(), k);
    writer.write(begin, end);
}

template<typename XS, typename YS, typename RX, typename RY, typename W>
void write_rays(const XS& xs, const YS& ys, const RX& rx, const RY& ry, const size_t& n, const size_t& k, W& writer) {
    writer.write(xs);
    write_strided(ys, k, writer);
    for (size_t i = 0; i < n; ++i)
        for (const auto& [x, y] : acstc::utils::zip(rx[i].data(), ry[i].data())) {
            write_strided(x, k, writer);
            write_strided(y, k, writer);
        }
}

template<typename KJ>
void save_rays(const acstc::config<types::real_t>& config, const std::string& filename, const bool binary, const KJ& k_j, const size_t& k) {
    const auto na = config.na();
    const auto nl = config.nl();
    const auto nm = k_j.size();

    const auto [rx, ry] = acstc::rays::compute(
        config.x0(), config.y_s(), config.l1(), nl, config.a0(), config.a1(), na, k_j);

    const auto& as = rx.template get<0>();
    const auto& ls = rx.template get<1>();
    const auto cl = nl / k;
    if (binary) {
        acstc::utils::binary_writer<types::real_t> writer(filename);
        writer.stream().write(reinterpret_cast<const char*>(&na), sizeof(na));
        writer.stream().write(reinterpret_cast<const char*>(&cl), sizeof(cl));
        writer.stream().write(reinterpret_cast<const char*>(&nm), sizeof(nm));
        write_rays(as, ls, rx, ry, nm, k, writer);
        return;
    }
    acstc::utils::text_writer<types::real_t> writer(filename);
    writer.stream() << as.size() << ' ' << cl << ' ' << nm << '\n';
    write_rays(as, ls, rx, ry, nm, k, writer);
}

template<typename V, typename W>
void write_conditions(const V& values, const size_t& size, W& writer) {
    for (const auto& it : values)
        writer.write(reinterpret_cast<const types::real_t*>(it.data()), 2 * size);
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
            ("config,c", po::value(&config_filename)->default_value("config.json"), "Config filename");

    po::options_description output("Output options");
    size_t step;
    std::string output_filename;
    output.add_options()
            ("output,o", po::value(&output_filename)->default_value("output.txt"), "Output filename")
            ("step,s", po::value(&step)->default_value(100)->value_name("k"), "Output every k-th computed row")
            ("binary", "Use binary output");

    po::options_description computation("Computation options");
    size_t num_workers, buff_size;
    computation.add_options()
            ("workers,w", po::value(&num_workers)->default_value(1), "Number of workers for computation")
            ("buff,b", po::value(&buff_size)->default_value(100), "Buff size to be used during multithreaded computation");

    po::options_description options;
    std::string job_type;
    options.add_options()
            ("job_type", po::value(&job_type)->default_value("solution"));
    options.add(generic).add(output).add(computation);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).positional(positional).options(options).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        po::options_description desc;
        desc.add(generic).add(output).add(computation);
        std::cout << "Usage: [solution|modes] (=solution) [options]\n" << desc << std::endl;
        return 0;
    }

    acstc::config config(config_filename);

    if (job_type == "solution") {
        acstc::solver solver(config);

        const auto [k0, phi_s] = config.create_source_modes();

        auto execute_function = [&](auto&& with_solver) {
            auto with_verbosity = add_verbosity(report, config.nx(), with_solver);
            auto with_writer = add_writer(config, output_filename, vm.count("binary") > 0, step, with_verbosity);
            with_writer();
        };

        auto with_solver = [&](const auto& init, const auto& k_j, const auto& phi_j, auto&& callback) mutable {
            solver.solve(init, k0, k_j, phi_j, callback, config.border_width(), config.past_n(), num_workers, buff_size);
        };

        auto with_solver_const = [&](const auto& init, const auto& k_j, const auto& phi_j, auto&& callback) mutable {
            solver.solve(init, k0, k_j, phi_j, callback, config.past_n(), num_workers, buff_size);
        };

        auto with_initial_conditions = add_initial_conditions(config, k0, phi_s, with_solver);
        auto with_initial_conditions_const = add_initial_conditions(config, k0, phi_s, with_solver_const);

        if (config.complex_modes())
            execute_function(add_modes<types::complex_t>(config, k0.size(), with_initial_conditions, with_initial_conditions_const));
        else
            execute_function(add_modes<types::real_t>(config, k0.size(), with_initial_conditions, with_initial_conditions_const));
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

    if (job_type == "rays") {
        if (config.const_modes()) {
            const auto [k_j, phi_j] = config.create_const_modes<types::real_t>(acstc::utils::verbosity::instance().level >= 2);
            save_rays(config, output_filename, vm.count("binary") > 0, k_j, step);
        } else {
            const auto [k_j, phi_j] = config.create_modes<types::real_t>(acstc::utils::verbosity::instance().level >= 2);
            save_rays(config, output_filename, vm.count("binary") > 0, k_j, step);
        }

        return 0;
    }

    if (job_type == "init") {
        const auto [k0, phi_s] = config.create_source_modes();
        const auto init = get_initial_conditions(config, k0, phi_s);

        const auto nj = k0.size();
        const auto ny = config.ny();

        const auto ys = acstc::utils::mesh_1d(config.y0(), config.y1(), ny);

        if (vm.count("binary") > 0) {
            acstc::utils::binary_writer<types::real_t> writer(output_filename);

            writer.stream().write(reinterpret_cast<const char*>(&nj), sizeof(nj));
            writer.stream().write(reinterpret_cast<const char*>(&ny), sizeof(ny));
            writer.write(ys);
            write_conditions(init, ny, writer);

            return 0;
        }

        acstc::utils::text_writer<types::real_t> writer(output_filename);

        writer.stream() << nj << ' ' << ny << '\n';
        writer.write(ys);
        write_conditions(init, ny, writer);

        return 0;
    }
    throw std::logic_error(std::string("Unknown job type: ") + job_type);
}
