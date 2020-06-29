#ifdef _WIN32
#define NOMINMAX
#endif

#include <cmath>
#include <chrono>
#include <complex>
#include <string>
#include <cstring>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "rays.hpp"
#include "config.hpp"
#include "solver.hpp"
#include "io/writer.hpp"
#include "utils/fft.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "nlohmann/json.hpp"
#include "utils/callback.hpp"
#include "utils/verbosity.hpp"
#include "initial_conditions.hpp"
#include "utils/progress_bar.hpp"
#include "utils/interpolation.hpp"
#include "boost/program_options.hpp"

namespace types = acstc::types;
namespace po = boost::program_options;

using nlohmann::json;
using acstc::utils::verboseln;
using acstc::utils::verboseln_lv;
using namespace std::complex_literals;

acstc::config<types::real_t> config;

enum class field_group : uint64_t {
    Modes = 1 << 0,
    Solver = 1 << 1,
    Rays = 1 << 2,
    Initial = 1 << 3
};

auto operator|(const field_group& a, const field_group& b) {
    return field_group(uint64_t(a) | uint64_t(b));
}

bool operator&(const field_group& a, const field_group& b) {
    return uint64_t(a) & uint64_t(b);
}

inline bool verbose(const size_t& level) {
    return acstc::utils::verbosity::instance().level >= level;
}

class to_string_helper {

public:

    explicit to_string_helper() {
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

    std::stringstream _stream;

} helper;

template<typename V>
auto get_strings(const V& values) {
    types::vector1d_t<std::string> result;

    for (const auto& it : values)
        result.push_back(helper.to_string(it));

    return result;
}

template<typename X, typename Y, typename V>
void print_table(const X& x, const Y& y, const V& values, std::stringstream& stream) {
    const auto sx = get_strings(x);
    const auto sy = get_strings(y);

    types::vector2d_t<std::string> sv;
    for (const auto& it : values)
        sv.emplace_back(get_strings(it));

    types::vector1d_t<int> widths(sy.size() + 1, 0);
    for (size_t i = 0; i < sy.size(); ++i) {
        widths[i + 1] = static_cast<int>(sy[i].size());
        for (size_t j = 0; j < sx.size(); ++j)
            widths[i + 1] = std::max(widths[i + 1], static_cast<int>(sv[j][i].size()));
    }

    for (const auto& it : sx)
        widths[0] = std::max(widths[0], static_cast<int>(it.size()));

    const size_t length = std::accumulate(widths.begin(), widths.end(), size_t(0)) + 3 * widths.size() + 1;

    char* sep_line = new char[length + 1];
    sep_line[length] = 0;

    std::memset(sep_line, '-', length);

    char* current = sep_line;
    for (const auto& it : widths) {
        *current = '+';
        current += it + 3;
    }
    *current = '+';

    stream << "        " << sep_line << '\n';

    char* buffer = new char[length + 1];
    current = buffer;
    current += sprintf(current, "| %*s |", widths[0], "");

    for (size_t i = 0; i < sy.size(); ++i)
        current += sprintf(current, " %*s |", widths[i + 1], sy[i].c_str());

    stream << "        " << buffer << '\n';

    *sep_line = '|';
    sep_line[length - 1] = '|';

    stream << "        " << sep_line << '\n';

    for (size_t i = 0; i < sv.size(); ++i) {
        current = buffer;
        current += sprintf(current, "| %*s |", widths[0], sx[i].c_str());
        for (size_t j = 0; j < sy.size(); ++j)
            current += sprintf(current, " %*s  ", widths[j + 1], sv[i][j].c_str());
        *--current = '|';
        stream << "        " << buffer << '\n';
    }

    *sep_line = '+';
    sep_line[length - 1] = '+';

    stream << "        " << sep_line << '\n';

    delete[] buffer;
    delete[] sep_line;
}

void print_loaded_data(const json& data, std::stringstream& stream) {
    const auto type = data[0].get<std::string>();

    if (type == "values") {
        stream << "loaded from config file;\n";
        return;
    }

    stream << "loaded from " << (type == "text_file" ? "text" : "binary") << " file " << data[1] << ";\n";
}

void print_bathymetry(std::stringstream& stream) {
    stream << "    Bathymetry: ";
    if (config.bathymetry().x().size() <= 5 && config.bathymetry().y().size() <= 5) {
        stream << '\n';
        print_table(config.bathymetry().x(), config.bathymetry().y(), config.bathymetry().data(), stream); 
        return;
    }

    print_loaded_data(config.data()["bathymetry"], stream);
}

template<typename T>
void print_array(const std::string& title, const json& data, std::stringstream& stream) {
    stream << "    " << title << ':';

    for (const auto& it : data)
        stream << " " << helper.to_string(it.get<T>());

    stream << ";\n";
}

void print_source_modes(const json& data, std::stringstream& stream) {
    if (data.count("k0") && data.count("phi_s")) {
        print_array<types::real_t>("Reference wave numbers", data["k0"], stream);
        print_array<types::real_t>("Source modal functions", data["phi_s"], stream);
    }
}

void print_mesh_spec(const std::string& title, const json& a, const json& b, const json& n, std::stringstream& stream) {
    stream << "    " << title << ": " << helper.to_string(a.get<types::real_t>()) << " <-- " << n << " --> " <<
        helper.to_string(b.get<types::real_t>()) << ";\n";
}

template<typename T>
void print_field(const std::string& title, const std::string& field, std::stringstream& stream) {
    stream << "    " << title << ": " << helper.to_string(config.data()[field].get<T>()) << ";\n";
}

void print_modes(std::stringstream& stream) {
    const auto& data = config.data();

    stream << "Modal parameters.\n";
    if (data.count("modes")) {
        stream << "    Modes: ";
        print_loaded_data(data["modes"], stream);
        print_source_modes(data, stream);
        stream << '\n';
        return;
    }

    print_field<double>("Mode subset", "mode_subset", stream);
    print_field<size_t>("Points per meter over z", "ppm", stream);
    print_field<size_t>("Richardson extrapolation order", "ordRich", stream);

    if (config.f_size() == 1)
        stream << "    Source frequency, Hz: " << config.f() << '\n';
    else {
        stream << "    Source function: ";
        if (config.f_size() <= 5)
            print_table(std::array<const char*, 1>{"f(t)"}, config.t(),
                std::array<std::remove_reference_t<decltype(config.source_function())>, 1>{config.source_function()}, stream);
        else
            print_loaded_data(data["source_function"], stream);
    }

    print_field<types::real_t>("Points per meter over z", "ppm", stream);
    print_field<types::real_t>("Source depth, m", "z_s", stream);
    print_field<bool>("Complex-valued modes", "complex_modes", stream);
    print_field<bool>("x-independent modes", "const_modes", stream);
    print_field<bool>("Use additive depth", "additive_depth", stream);

    stream << "    Bottom layers (top_speed -> bottom_speed; depth; density):\n";
    for (const auto& [c1, c2, z, r] : acstc::utils::zip(
        data["bottom_c1s"],
        data["bottom_c2s"],
        data["bottom_layers"],
        data["bottom_rhos"]
    ))
        stream << "        " << helper.to_string(c1.get<types::real_t>()) << " -> " << 
                                helper.to_string(c2.get<types::real_t>()) << "; " << 
                                helper.to_string( z.get<types::real_t>()) << "; " << 
                                helper.to_string( r.get<types::real_t>()) << ";\n";

    stream << "    Number of water layers: " << data["n_layers"] << ";\n";
    print_array<types::real_t>("Beta parameters", data["betas"], stream);

    const auto nm = data["max_mode"].get<size_t>();
    stream << "    Maximal number of modes: " << (nm == size_t(-1) ? "All" : helper.to_string(nm)) << ";\n";

    const auto mn = data["n_modes"].get<size_t>();
    stream << "    Required number of modes: " << (mn == size_t(-1) ? "All" : helper.to_string(mn)) << ";\n";

    print_bathymetry(stream);

    stream << "    Hydrology: ";
    print_loaded_data(data["hydrology"], stream);

    if (config.const_modes())
        stream << "    Number of points over y: " << 
            (data.count("mny") ? data["mny"].get<size_t>() : config.bathymetry().y().size()) << ";\n";
    else {
        size_t nx, ny;
        if (data.count("mnx") && data.count("mny")) {
            nx = data["mnx"].get<size_t>();
            ny = data["mny"].get<size_t>();
        } else {
            nx = config.bathymetry().x().size();
            ny = config.bathymetry().y().size();
        }
        stream << "    Number of points over x: " << nx << ";\n    Number of points over y: " << ny << ";\n";
    }

    print_source_modes(data, stream);

    stream << '\n';
}

void print_tapering(std::stringstream& stream) {
    const auto& [type, vl, vr] = config.get_tapering_parameters();
    stream << "    Tapering:\n        Type: " << type << ";\n" << 
        "        Left: " << vl << ";\n" <<
        "        Right: " << vr << ";\n";
}

void print_initial_conditions(std::stringstream& stream) {
    const auto& data = config.data();

    stream << "Initial conditions parameters.\n    Type: ";
    const auto type = data["init"].get<std::string>();

    if (type == "greene")
        stream << "Greene source";
    else if (type == "gauss")
        stream << "Gaussian source";
    else if (type == "ray_simple") {
        stream << "Ray-based source assuming homogeneous medium;\n";
        print_mesh_spec("Angle mesh", data["a0"], data["a1"], data["na"], stream);
        print_tapering(stream);
        stream << '\n';
        return;
    }
    else if (type == "ray") {
        stream << "Ray-based source;\n";
        print_mesh_spec("Angle mesh", data["a0"], data["a1"], data["na"], stream);
        print_mesh_spec("Natural parameter mesh", data["l0"], data["l1"], data["nl"], stream);
        print_tapering(stream);
        stream << '\n';
        return;
    }
    else
        throw std::logic_error(std::string("Unknown initial conditions type: ") + type);

    stream << ";\n\n";
}

void print_solver(std::stringstream& stream) {
    const auto& data = config.data();

    stream << "Solver parameters:\n";

    const auto& receivers = config.receiver_depth();
    stream << "    Receivers:\n";
    for (size_t i = 0; i < std::min(size_t(5), receivers.data().size()); ++i) {
        const auto& [x, y] = receivers.points()[i];
        stream << "        [ " << x << ", " << y << ", " << receivers.data()[i] << " ];\n";
    }

    print_field<types::real_t>("Source z coordinate", "y_s", stream);
    print_field<size_t>("Width of smooth border over edges", "border_width", stream);

    const auto pn = data["past_n"].get<size_t>();
    stream << "    Size of border convolution: "  << (pn == 0 ? "Full" : helper.to_string(pn)) << ";\n";

    print_mesh_spec("x mesh", data["x0"], data["x1"], data["nx"], stream);
    print_mesh_spec("y mesh", data["y0"], data["y1"], data["ny"], stream);

    stream << "    Root approximation coeffitients:\n" << 
        "        a: " << helper.to_string(config.a()) << ";\n" <<
        "        b: " << helper.to_string(config.b()) << ";\n" <<
        "        c: " << helper.to_string(config.c()) << ";\n";

    stream << '\n';
}

void verbose_config_field_group_parameters(const field_group& group) {
    if (!verbose(3))
        return;

    std::stringstream stream;

    if (group & field_group::Modes)
        print_modes(stream);

    if (group & field_group::Solver)
        print_solver(stream);

    if (group & field_group::Initial)
        print_initial_conditions(stream);

    if (group & field_group::Rays) {
        const auto& data = config.data();

        stream << "Rays parameters:\n";
        print_mesh_spec("Angle mesh", data["a0"], data["a1"], data["na"], stream);
        print_mesh_spec("Natural parameter mesh", data["l0"], data["l1"], data["nl"], stream);
    }

    std::cout << stream.str();
}

template<typename F>
auto pass_tapering(const F& func) {
    const auto& [type, vl, vr] = config.get_tapering_parameters();
    if (type == "percentage")
        return func(acstc::percentage_tapering(vl, vr));
    if (type == "angled")
        return func(acstc::angled_tapering(vl, vr));

    throw std::logic_error(std::string("Unknown tapering type: ") + type);
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const KS& k0, const PS& phi_s,
    const acstc::utils::linear_interpolated_data_1d<types::real_t>& k_j) {
    return pass_tapering(
        [&](const auto& tapering) {
            return acstc::ray_source(config.x0(), config.y0(), config.y1(), config.ny(), 0., config.y_s(), config.l1(), config.nl(),
                config.a0(), config.a1(), config.na(), k0, phi_s, k_j, tapering);
        }
    );
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const KS& k0, const PS& phi_s) {
    auto [k_j, phi_j] = config.create_const_modes<types::real_t>(config.n_modes(), verbose(2));

    if (k_j.size() > k0.size())
        k_j.erase_last(k_j.size() - k0.size());

    return get_ray_initial_conditions(k0, phi_s, k_j);
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const KS& k0, const PS& phi_s,
    const acstc::utils::linear_interpolated_data_1d<types::real_t, types::complex_t>& k_j) {
    const auto& ys = k_j.get<0>();
    types::vector2d_t<types::real_t> new_k_j(k_j.size(), types::vector1d_t<types::real_t>(ys.size()));

    for (size_t j = 0; j < k_j.size(); ++j)
        std::transform(k_j[j].data().begin(), k_j[j].data().end(), new_k_j[j].begin(), [](const auto& v) { return v.real(); });

    return get_ray_initial_conditions(k0, phi_s, acstc::utils::linear_interpolated_data_1d<types::real_t>(ys, new_k_j));
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const KS& k0, const PS& phi_s,
    const acstc::utils::linear_interpolated_data_2d<types::real_t>& k_j) {
    types::vector2d_t<types::real_t> new_k_j;
    new_k_j.reserve(k_j.size());

    for (size_t j = 0; j < k_j.size(); ++j)
        new_k_j.emplace_back(k_j[j][0].begin(), k_j[j][0].end());

    return get_ray_initial_conditions(k0, phi_s, acstc::utils::linear_interpolated_data_1d<types::real_t>(k_j.get<1>(), new_k_j));
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const KS& k0, const PS& phi_s,
    const acstc::utils::linear_interpolated_data_2d<types::real_t, types::complex_t>& k_j) {
    const auto& ys = k_j.get<1>();
    types::vector2d_t<types::real_t> new_k_j(k_j.size(), types::vector1d_t<types::real_t>(ys.size()));

    for (size_t j = 0; j < k_j.size(); ++j)
        std::transform(k_j[j][0].begin(), k_j[j][0].end(), new_k_j[j].begin(), [](const auto& v) { return v.real(); });

    return get_ray_initial_conditions(k0, phi_s, acstc::utils::linear_interpolated_data_1d<types::real_t>(ys, new_k_j));
}

template<typename KS, typename PS>
auto get_simple_initial_conditions(const KS& k0, const PS& phi_s) {
    const auto init = config.init();

    KS ws(k0.size());
    PS as(phi_s.size());
    std::transform(phi_s.begin(), phi_s.end(), as.begin(), [](const auto& phi) { return phi / (2 * std::sqrt(M_PI)); });
    std::transform(k0.begin(), k0.end(), ws.begin(), [](const auto& k0) { return 1 / std::pow(k0, 2); } );

    if (init == "greene")
        return acstc::greene_source<types::complex_t>(config.y0(), config.y1(), config.ny(), config.y_s(), as, ws);

    if (init == "gauss")
        return acstc::gaussian_source<types::complex_t>(config.y0(), config.y1(), config.ny(), config.y_s(), as, ws);

    if (init == "ray_simple")
        return pass_tapering(
            [&](const auto& tapering) {
                return acstc::simple_ray_source(config.x0(), config.y0(), config.y1(), config.ny(), 
                    config.a0(), config.a1(), k0, phi_s, tapering);
            }
        );

    throw std::logic_error(std::string("Unknown initial conditions type: ") + init);
}

template<typename KS, typename PS, typename KJ>
auto get_initial_conditions(const KS& k0, const PS& phi_s, const KJ& k_j) {
    const auto init = config.init();

    if (init == "ray")
        return get_ray_initial_conditions(k0, phi_s, k_j);

    return get_simple_initial_conditions(k0, phi_s);
}

template<typename KS, typename PS>
auto get_initial_conditions(const KS& k0, const PS& phi_s) {
    const auto init = config.init();

    if (init == "ray")
        return get_ray_initial_conditions(k0, phi_s);

    return get_simple_initial_conditions(k0, phi_s);
}

template<typename F>
auto add_verbosity(const size_t& report, const size_t& n, F& function) {
    return [&](auto&& writer) mutable {
        if (verbose(2))
            function(acstc::utils::callbacks(writer, acstc::utils::progress_bar_callback(n, "Solution", true)));
        else {
            const auto start = std::chrono::system_clock::now();
            function(writer);
            const auto end = std::chrono::system_clock::now();
            verboseln(1, "Elapsed time: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");
        }
        return;
    };
}

template<typename F>
auto add_writer(const std::string& filename,
        const bool binary, const size_t& step, F& function) {
    return [&, binary]() mutable {
        const auto xs = acstc::utils::mesh_1d(config.x0(), config.x1(), (config.nx() - 1) / step + 1);
        const auto ys = acstc::utils::mesh_1d(config.y0(), config.y1(), config.ny());
        if (binary) {
            const auto nf = static_cast<const uint32_t>(config.f_size());
            const auto nx = static_cast<const uint32_t>(xs.size());
            const auto ny = static_cast<const uint32_t>(ys.size());
            acstc::utils::binary_writer<types::real_t> writer(filename);
            writer.stream().write(reinterpret_cast<const char*>(&nf), sizeof(nf));
            writer.stream().write(reinterpret_cast<const char*>(&nx), sizeof(nx));
            writer.stream().write(reinterpret_cast<const char*>(&ny), sizeof(ny));
            writer.write(xs);
            writer.write(ys);
            function(acstc::utils::ekc_callback(step, [&](const auto& x, const auto& data) {
                writer.write(reinterpret_cast<const types::real_t*>(data.data()), ny * 2);
            }));
            return;
        }
        acstc::utils::text_writer<types::real_t> writer(filename);
        writer.stream() << config.f_size() << ' ' << xs.size() << ' ' << ys.size() << '\n';
        writer.write(xs);
        writer.write(ys);
        function(acstc::utils::ekc_callback(step, [&](const auto& x, const auto& data) {
            writer.write(reinterpret_cast<const types::real_t*>(data.data()), ys.size() * 2);
        }));
    };
}

template<typename T, typename F1, typename F2>
auto add_modes(const size_t mn, F1& function, F2& function_const) {
    return [&, mn](auto&& callback) mutable {
        if (config.const_modes()) {
            const auto start = std::chrono::system_clock::now();
            auto [k_j, phi_j] = config.create_const_modes<T>(mn, verbose(2));
            const auto end = std::chrono::system_clock::now();
            verboseln_lv(1, "Modes computing time: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");

            if (k_j.size() > mn) {
                k_j.erase_last(k_j.size() - mn);
                phi_j.erase_last(phi_j.size() - mn);
            }
            function_const(k_j, phi_j, callback);
            return;
        }

        const auto start = std::chrono::system_clock::now();
        auto [k_j, phi_j] = config.create_modes<T>(mn, verbose(2));
        const auto end = std::chrono::system_clock::now();
        verboseln_lv(1, "Modes computing time: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");

        if (k_j.size() > mn) {
            k_j.erase_last(k_j.size() - mn);
            phi_j.erase_last(phi_j.size() - mn);
        }
        function(k_j, phi_j, callback);
    };
}

template<typename KS, typename PS, typename F>
auto add_initial_conditions(const KS& k0, const PS& phi_s, F& function) {
    return [&](const auto& k_j, const auto& phi_j, auto&& callback) mutable {
        function(get_initial_conditions(k0, phi_s, k_j), k_j, phi_j, callback);
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
void save_modes(const std::string& filename, const bool binary, const W& wrapper) {
    if (config.const_modes()) {
        const auto [k_j, phi_j] = config.create_const_modes<T>(config.n_modes(), acstc::utils::verbosity::instance().level >= 2);
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
    const auto [k_j, phi_j] = config.create_modes<T>(config.n_modes(), acstc::utils::verbosity::instance().level >= 2);
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
void save_rays(const std::string& filename, const bool binary, const KJ& k_j, const size_t& k) {
    const auto na = static_cast<uint32_t>(config.na());
    const auto nl = static_cast<uint32_t>(config.nl());
    const auto nm = static_cast<uint32_t>(k_j.size());

    const auto [rx, ry] = acstc::rays::compute(
        config.x0(), config.y_s(), config.l1(), nl, config.a0(), config.a1(), na, k_j, verbose(2));

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
    writer.stream() << na << ' ' << cl << ' ' << nm << '\n';
    write_rays(as, ls, rx, ry, nm, k, writer);
}

template<typename V, typename W>
void write_conditions(const V& values, const size_t& size, W& writer) {
    for (const auto& it : values)
        writer.write(reinterpret_cast<const types::real_t*>(it.data()), 2 * size);
}

template<typename W>
void write_impulse(const types::vector1d_t<types::real_t>& tau,
    const types::vector2d_t<types::real_t>& impulse, W&& writer) {
    writer.write(config.t());
    writer.write(tau);
    for (const auto& it : impulse)
        writer.write(it);
}

void save_impulse(const std::string& filename, const types::vector1d_t<types::real_t>& tau,
                  const types::vector2d_t<types::real_t>& impulse, const bool& binary) {
    const auto nr = static_cast<uint32_t>(tau.size());
    const auto nt = static_cast<uint32_t>(impulse[0].size());
    if (binary) {
        acstc::utils::binary_writer<types::real_t> writer(filename);
        writer.stream().write(reinterpret_cast<const char*>(&nr), sizeof(nr));
        writer.stream().write(reinterpret_cast<const char*>(&nt), sizeof(nt));
        write_impulse(tau, impulse, writer);
        return;
    }
    acstc::utils::text_writer<types::real_t> writer(filename);
    writer.stream() << nr << ' ' << nt << '\n';
    write_impulse(tau, impulse, writer);
}

int main(int argc, char* argv[]) {
    try {
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

        config.update_from_file(config_filename);

        if (job_type == "solution") {
            verbose_config_field_group_parameters(field_group::Modes | field_group::Solver | field_group::Initial);
            acstc::utils::progress_bar::leave = true;

            const auto nf = config.f_size();
            acstc::solver solver(config);
            types::vector1d_t<types::real_t> k0, phi_s;

            auto with_solver = [&](const auto& init, const auto& k_j, const auto& phi_j, auto&& callback) mutable {
                solver.solve(init, k0, k_j, phi_j, callback, config.border_width(), config.past_n(), num_workers, buff_size);
            };

            auto with_solver_const = [&](const auto& init, const auto& k_j, const auto& phi_j, auto&& callback) mutable {
                solver.solve(init, k0, k_j, phi_j, callback, config.past_n(), num_workers, buff_size);
            };

            auto with_initial_conditions = add_initial_conditions(k0, phi_s, with_solver);
            auto with_initial_conditions_const = add_initial_conditions(k0, phi_s, with_solver_const);

            auto function = [&](auto&& writer) {
                acstc::utils::progress_bar pbar(nf, "Frequency", verbose(2), acstc::utils::progress_bar::leave_behavior::stay);

                auto execute_function = [&](auto&& with_solver) {
                    add_verbosity(report, config.nx(), with_solver)(writer);
                };

                for (const auto& fi : pbar) {
                    config.f_index(fi);
                    std::tie(k0, phi_s) = config.create_source_modes(config.n_modes());

                    if (config.complex_modes())
                        execute_function(add_modes<types::complex_t>(k0.size(), with_initial_conditions, with_initial_conditions_const));
                    else
                        execute_function(add_modes<types::real_t>(k0.size(), with_initial_conditions, with_initial_conditions_const));
                }

            };

            add_writer(output_filename, vm.count("binary") > 0, step, function)();
            return 0;
        }

        if (job_type == "impulse") {
            verbose_config_field_group_parameters(field_group::Modes | field_group::Solver | field_group::Initial);
            acstc::utils::progress_bar::leave = true;

            acstc::utils::real_fft<types::real_t, types::complex_t> fft(static_cast<int>(config.f_size()));
            std::memset(fft.backward_data(), 0, sizeof(types::complex_t) * fft.size());
            std::memcpy(fft.forward_data(), config.source_function().data(), fft.size() * sizeof(types::real_t));
            fft.execute_forward();
            std::transform(std::as_const(fft).backward_data(), fft.backward_data_end(), fft.backward_data(), [](const auto& v) { return std::conj(v); });

            config.f_mode(decltype(config)::mode::impulse);
            const auto nf = config.f_size();

            acstc::solver solver(config);

            const auto& depth = config.receiver_depth();
            const auto nr = depth.points().size();

            types::vector1d_t<types::real_t> k0, phi_s;
            types::vector2d_t<types::complex_t> result(nr, types::vector1d_t<types::complex_t>(nf, types::real_t(0)));

            types::vector1d_t<size_t> ix(nr);
            std::iota(ix.begin(), ix.end(), 0);
            std::sort(ix.begin(), ix.end(), [&](const auto& a, const auto& b) {
                return std::get<0>(depth.points(a)) < std::get<0>(depth.points(b));
                });

            const auto iy = acstc::utils::mesh_1d(config.y0(), config.y1(), config.ny());

            size_t li = 0;
            acstc::utils::progress_bar pbar(nf, "Impulse", verbose(2), acstc::utils::progress_bar::leave_behavior::stay);

            auto execute_function = [&, &fi = pbar.current()](auto&& with_solver) {
                with_solver(acstc::utils::callbacks(
                    acstc::utils::progress_bar_callback(config.nx(), "Solution", verbose(2)),
                    [&, last = types::vector1d_t<types::complex_t>(), last_x = double(0)](const auto& mx, const auto& data) mutable {
                    const auto x = mx + config.x0();
                    for (; last.size() && li < ix.size() && std::get<0>(depth.points(ix[li])) <= x; ++li) {
                        const auto& [px, py] = depth.points(ix[li]);
                        const auto& [ya, yb] = acstc::utils::find_indices(iy, py);
                        result[ix[li]][fi] = acstc::utils::__impl::linear_interpolation::field_point(
                            last[ya], last[yb], data[ya], data[yb], last_x, x, iy[ya], iy[yb], px, py) * fft.backward_data(fi);
                    }
                    last_x = x;
                    last.assign(data.begin(), data.end());
                }));
            };

            auto with_solver = [&](const auto& init, const auto& k_j, const auto& phi_j, auto&& callback) mutable {
                solver.solve(init, k0, k_j, phi_j, callback, config.border_width(), config.past_n(), num_workers, buff_size);
            };

            auto with_solver_const = [&](const auto& init, const auto& k_j, const auto& phi_j, auto&& callback) mutable {
                solver.solve(init, k0, k_j, phi_j, callback, config.past_n(), num_workers, buff_size);
            };

            auto with_initial_conditions = add_initial_conditions(k0, phi_s, with_solver);
            auto with_initial_conditions_const = add_initial_conditions(k0, phi_s, with_solver_const);

            const auto max = std::abs(*std::max_element(std::as_const(fft).backward_data(), fft.backward_data_end(),
                [](const auto& a, const auto& b) {
                    return std::abs(a) < std::abs(b);
                }
            ));

            for (const auto& fi : pbar) {
                li = 0;
                config.f_index(fi);

                if (std::abs(fft.backward_data(fi)) < max * config.tolerance())
                    continue;

                std::tie(k0, phi_s) = config.create_source_modes(0);

                if (config.complex_modes())
                    execute_function(add_modes<types::complex_t>(k0.size(), with_initial_conditions, with_initial_conditions_const));
                else
                    execute_function(add_modes<types::real_t>(k0.size(), with_initial_conditions, with_initial_conditions_const));
            }

            const auto h = config.bathymetry().point(0, config.y_s());
            const auto cs = config.hydrology().line(0, 0, h, static_cast<size_t>(h));
            const auto cm = std::max({
                *std::max_element(cs.begin(), cs.end()),
                *std::max_element(config.bottom_c1s().begin(), config.bottom_c1s().end()),
                *std::max_element(config.bottom_c2s().begin(), config.bottom_c2s().end())
                });

            auto omeg = acstc::utils::mesh_1d(0., 1 / config.dt(), fft.size());
            std::transform(omeg.begin(), omeg.end(), omeg.begin(), [](const auto& value) { return value * 2 * M_PI; });

            types::vector1d_t<types::real_t> tau(nr);
            types::vector2d_t<types::real_t> impulse(nr, types::vector1d_t<types::real_t>(fft.size()));

            for (size_t i = 0; i < result.size(); ++i) {
                const auto d = tau[i] = std::hypot(std::get<0>(depth.points(i)), std::get<1>(depth.points(i))) / cm;

                auto zip = acstc::utils::zip(result[i], omeg);
                std::transform(zip.begin(), zip.end(), fft.backward_data(),
                    [&](const auto& v) {
                        return std::conj(std::get<0>(v) * std::exp(-1i * d * std::get<1>(v)));
                    }
                );

                fft.execute_backward().normalize_forward();
                std::memcpy(impulse[i].data(), fft.forward_data(), fft.size() * sizeof(*fft.forward_data()));
            }

            save_impulse(output_filename, tau, impulse, vm.count("binary") > 0);

            return 0;
        }

        if (job_type == "modes") {
            verbose_config_field_group_parameters(field_group::Modes);

            if (config.complex_modes())
                save_modes<types::complex_t>(output_filename, vm.count("binary") > 0,
                    [](auto& writer) {
                        return [&](const auto& data) {
                            writer.write(reinterpret_cast<const types::real_t*>(data.data()), data.size() * 2);
                        };
                    });
            else
                save_modes<types::real_t>(output_filename, vm.count("binary") > 0,
                    [](auto& writer) { return [&](const auto& data) { writer(data); }; });
            return 0;
        }

        if (job_type == "rays") {
            verbose_config_field_group_parameters(field_group::Rays);

            if (config.const_modes()) {
                const auto [k_j, phi_j] = config.create_const_modes<types::real_t>(config.n_modes(), verbose(2));
                save_rays(output_filename, vm.count("binary") > 0, k_j, step);
            }
            else {
                const auto [k_j, phi_j] = config.create_modes<types::real_t>(config.n_modes(), verbose(2));
                save_rays(output_filename, vm.count("binary") > 0, k_j, step);
            }

            return 0;
        }

        if (job_type == "init") {
            verbose_config_field_group_parameters(field_group::Initial);

            const auto [k0, phi_s] = config.create_source_modes(config.n_modes());
            const auto init = get_initial_conditions(k0, phi_s);

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
    catch (const std::exception& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
}
