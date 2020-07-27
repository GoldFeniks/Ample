#ifdef _WIN32
#define NOMINMAX
#endif

#include <set>
#include <cmath>
#include <chrono>
#include <string>
#include <complex>
#include <cstring>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <filesystem>
#include <functional>
#include <unordered_set>
#include "rays.hpp"
#include "config.hpp"
#include "solver.hpp"
#include "io/writer.hpp"
#include "utils/fft.hpp"
#include "feniks/zip.hpp"
#include "utils/types.hpp"
#include "utils/utils.hpp"
#include "nlohmann/json.hpp"
#include "utils/callback.hpp"
#include "utils/verbosity.hpp"
#include "boost/lexical_cast.hpp"
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
    Nothing = 0,
    Modes = 1 << 1,
    Solver = 1 << 2,
    Rays = 1 << 3,
    Initial = 1 << 4
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

template<typename V, typename W>
void write_modes(const acstc::utils::linear_interpolated_data_1d<types::real_t, V>& modes, W&& writer) {
    for (size_t i = 0; i < modes.size(); ++i)
        writer(modes[i].data());
}

template<typename V, typename W>
void write_modes(const acstc::utils::linear_interpolated_data_2d<types::real_t, V>& modes, W&& writer) {
    for (size_t i = 0; i < modes.size(); ++i)
        for (const auto& it : modes[i].data())
            writer(it);
}

template<typename KV, typename PV, typename W>
void write_modes(const KV& k_j, const PV& phi_j, W&& writer) {
    write_modes(k_j, writer);
    write_modes(phi_j, writer);
}

template<typename W>
auto wrap_writer(W& writer) {
    return [&writer](const auto& data) {
        writer.write(reinterpret_cast<const types::real_t*>(data.data()), data.size() * sizeof(data[0]) / sizeof(types::real_t));
    };
}

template<typename V, typename W>
void write_strided(const V& v, const size_t& k, W&& writer) {
    auto [begin, end] = acstc::utils::stride(v.begin(), v.end(), k);
    writer.write(begin, end);
}

template<typename RX, typename RY, typename W>
void write_rays(const RX& rx, const RY& ry, const size_t& n, const size_t& k, W&& writer) {
    for (size_t i = 0; i < n; ++i)
        for (const auto& [x, y] : feniks::zip(rx[i].data(), ry[i].data())) {
            write_strided(x, k, writer);
            write_strided(y, k, writer);
        }
}

template<typename KJ>
void save_rays(const std::string& filename, const bool binary, const KJ& k_j, const size_t& k) {
    const auto na = config.na();
    const auto nl = config.nl();
    const auto nm = k_j.size();

    const auto [rx, ry] = acstc::rays::compute(
        config.x0(), config.y_s(), config.l1(), nl, config.a0(), config.a1(), na, k_j, verbose(2));

    if (binary)
        write_rays(rx, ry, nm, k, acstc::utils::binary_writer<types::real_t>(filename));
    else
        write_rays(rx, ry, nm, k, acstc::utils::text_writer<types::real_t>(filename));
}

template<typename V, typename W>
void write_conditions(const V& values, W&& writer) {
    for (const auto& it : values)
        writer.write(reinterpret_cast<const types::real_t*>(it.data()), 2 * it.size());
}

template<typename W>
void write_impulse(const types::vector2d_t<types::real_t>& impulse, W&& writer) {
    for (const auto& it : impulse)
        writer.write(it);
}

void save_impulse(const std::string& filename,
                  const types::vector2d_t<types::real_t>& impulse, const bool& binary) {
    if (binary)
        write_impulse(impulse, acstc::utils::binary_writer<types::real_t>(filename));
    else
        write_impulse(impulse, acstc::utils::text_writer<types::real_t>(filename));
}

template<typename S, typename K, typename V, typename I, typename C>
void solve(S& solver, const I& init, const K& k0,
           const acstc::utils::linear_interpolated_data_1d<types::real_t, V>& k_j, 
           const acstc::utils::linear_interpolated_data_1d<types::real_t>& phi_j,
           C&& callback, const size_t& num_workers, const size_t& buff_size) {
    solver.solve(init, k0, k_j, phi_j, callback, config.past_n(), num_workers, buff_size);
}

template<typename S, typename K, typename V, typename I, typename C>
void solve(S& solver, const I& init, const K& k0,
           const acstc::utils::linear_interpolated_data_2d<types::real_t, V>& k_j, 
           const acstc::utils::linear_interpolated_data_2d<types::real_t>& phi_j,
           C&& callback, const size_t& num_workers, const size_t& buff_size) {
    solver.solve(init, k0, k_j, phi_j, callback, config.border_width(), config.past_n(), num_workers, buff_size);
}

const std::set<std::string> available_jobs{
    "sel",
    "init",
    "rays",
    "modes",
    "impulse",
    "solution"
};

class jobs {

public:

    jobs() = default;
    jobs(std::initializer_list<std::string> init) : _jobs(init) {}

    void add_job(const std::string& job) {
        if (available_jobs.find(job) == available_jobs.end())
            throw std::logic_error(std::string("Unknown job type: ") + job);
        _jobs.insert(job);
    }

    bool has_job(const std::string& name) const {
        return _jobs.find(name) != _jobs.end();
    }

    auto size() const {
        return _jobs.size();
    }

    const auto& raw() const {
        return _jobs;
    }

private:

    std::unordered_set<std::string> _jobs;

};

class jobs_config {

public:

    ::jobs jobs;
    bool binary;
    size_t step, num_workers, buff_size;
    std::filesystem::path output, config_path;

    void perform() {
        config.update_from_file(config_path);
        acstc::utils::progress_bar::clear_on_end = true;
        std::filesystem::create_directories(output);

        field_group group = field_group::Nothing;

        _prep("init", field_group::Initial, group);
        _prep("rays", field_group::Rays, group);
        _prep("modes", field_group::Modes, group);
        _prep("solution", field_group::Modes | field_group::Solver | field_group::Initial, group);

        if (jobs.has_job("impulse"))
            group = group | field_group::Modes | field_group::Solver | field_group::Initial;

        verbose_config_field_group_parameters(group);

        _meta["step"] = step;
        _meta["binary"] = binary;
        _meta["output"] = output;
        _meta["f"] = json::array();
        _meta["k0"] = json::array();
        _meta["mnx"] = config.mnx();
        _meta["mny"] = config.mny();
        _meta["jobs"] = jobs.raw();
        _meta["phi_s"] = json::array();
        _meta["n_modes"] = json::array();
        _meta["num_workers"] = num_workers;
        _meta["original_config_path"] = config_path;
        _meta["n_receivers"] = config.receiver_depth().points().size();

        _add_mesh("x", config.x0(), config.x1(), (config.nx() - 1) / step + 1);
        _add_mesh("y", config.y0(), config.y1(), config.ny());
        _add_mesh("a", config.a0(), config.a1(), config.na());
        _add_mesh("l", config.l0(), config.l1(), (config.nl() - 1) / step + 1);

        _meta["t"] = {
            { "a", config.t().front() },
            { "b", config.t().back()  },
            { "n", config.t().size()  },
            { "d", config.dt() }
        };

        _pick_writer();

        std::ofstream out(output / "meta.json");
        out << std::setw(4) << _meta << std::endl;

        config.save(output);
    }

private:

    json _meta;

    void _prep(const char* name, const field_group& params, field_group& group) {
        if (jobs.has_job(name)) {
            std::filesystem::create_directories(output / name);
            group = group | params;
        }
    }

    auto _add_extension(std::filesystem::path path) const {
        path += binary ? ".bin" : ".txt";
        return path;
    }

    auto _get_filename(const char* name) const {
        return _add_extension(output / name / boost::lexical_cast<std::string>(config.f()));
    }

    void _add_mesh(const char* name, const types::real_t& a, const types::real_t& b, const size_t& n) {
        _meta[name] = {
            { "a", a },
            { "b", b },
            { "n", n },
            { "d", (b - a) / (n - 1)}
        };
    }

    void _pick_writer() {
        if (binary)
            _pick_const<acstc::utils::binary_writer<types::real_t>>();
        else
            _pick_const<acstc::utils::text_writer<types::real_t>>();
    }

    template<typename W>
    void _pick_const() {
        if (config.const_modes())
            _pick_complex<W, true>();
        else
            _pick_complex<W, false>();
    }

    template<typename W, bool Const>
    void _pick_complex() {
        if (config.complex_modes())
            performer<W, make_modes<Const, types::complex_t>>(*this).perform();
        else
            performer<W, make_modes<Const, types::real_t>>(*this).perform();

    }

    template<bool Const, typename T>
    struct make_modes {

        static auto make(const size_t& nm, const bool& show_progress) {
            if constexpr (Const)
                return config.create_const_modes<T>(nm, show_progress);
            else
                return config.create_modes<T>(nm, show_progress);
        }

    };

    template<typename W, typename M>
    class performer {

    public:

        performer(jobs_config& owner) : _owner(owner) {}

        void perform() {
            const auto has_sel = _owner.jobs.has_job("sel");
            const auto has_impulse = _owner.jobs.has_job("impulse");

            if (has_impulse) {
                 _fft = new acstc::utils::real_fft<types::real_t, types::complex_t>(static_cast<int>(config.f_size()));

                std::memset(_fft->backward_data(), 0, sizeof(types::complex_t) * _fft->size());
                std::memcpy(_fft->forward_data(), config.source_function().data(), _fft->size() * sizeof(types::real_t));
                _fft->execute_forward();

                std::transform(std::as_const(*_fft).backward_data(), _fft->backward_data_end(), _fft->backward_data(), 
                    [](const auto& v) { 
                        return std::conj(v); 
                    }
                );

                config.f_mode(decltype(config)::mode::impulse);

                const auto& depth = config.receiver_depth();
                const auto nr = depth.points().size();

                _impulse_result = new types::vector2d_t<types::complex_t>(nr,
                    types::vector1d_t<types::complex_t>(config.f_size(), types::real_t(0))
                );

                _ix = new types::vector1d_t<size_t>(nr);
                std::iota(_ix->begin(), _ix->end(), 0);
                std::sort(_ix->begin(), _ix->end(),
                    [&depth](const auto& a, const auto& b) {
                        return std::get<0>(depth.points(a)) < std::get<0>(depth.points(b));
                    }
                );

                _iy = new types::vector1d_t<types::real_t>(acstc::utils::mesh_1d(config.y0(), config.y1(), config.ny()));

                _max = std::abs(*std::max_element(std::as_const(*_fft).backward_data(), _fft->backward_data_end(),
                    [](const auto& a, const auto& b) {
                        return std::abs(a) < std::abs(b);
                    }
                ));
            }

            if (has_sel)
                _sel_result = new types::vector2d_t<types::real_t>(
                    config.nx() / _owner.step + 1,
                    types::vector1d_t<types::real_t>(config.ny(), 0)
                );

            const auto [f0, f1] = config.sel_range();
            acstc::utils::progress_bar pbar(config.f_size(), "Frequency", verbose(2), acstc::utils::progress_bar::on_end::leave);
            for (const auto& fi : pbar) {
                config.f_index(fi);
                pbar.set_description(std::string("Frequency ") + boost::lexical_cast<std::string>(config.f()));

                const auto f = config.f();
                if (has_impulse && std::abs(_fft->backward_data(fi)) < _max * config.tolerance() &&
                    !(has_sel && f0 <= f && f1 >= f))
                    continue;

                const auto [k0, phi_s] = config.create_source_modes(config.n_modes());
                if (!k0.size())
                    continue;

                _owner._meta["f"].push_back(config.f());
                const auto [k_j, phi_j] = _perform_modes(k0, phi_s);
                const auto init = _perform_init(k0, phi_s, k_j);

                _perform_rays(k_j, phi_j);

                _perform_sel(init, k0, k_j, phi_j);

                if (has_impulse) {
                    auto& result = *_impulse_result;
                    for (size_t i = 0; i < result.size(); ++i)
                        if (i != config.reference_index())
                            result[i][fi] *= _s;
                }
            }

            if (has_sel) {
                const auto size = _fft != nullptr ? _fft->size() : config.f_size();

                for (auto& x : *_sel_result)
                    for (auto& y : x)
                        y *= config.dt() / size;

                W writer(_owner._add_extension(_owner.output / "sel"));
                for (const auto& it : *_sel_result)
                    writer.write(it);

                delete _sel_result;
            }

            if (has_impulse) {
                const auto h = config.bathymetry().point(0, config.y_s());
                const auto cs = config.hydrology().line(0, 0, h, static_cast<size_t>(h));
                const auto cm = std::max({
                    *std::max_element(cs.begin(), cs.end()),
                    *std::max_element(config.bottom_c1s().begin(), config.bottom_c1s().end()),
                    *std::max_element(config.bottom_c2s().begin(), config.bottom_c2s().end())
                });

                auto omeg = acstc::utils::mesh_1d(0., 1 / config.dt(), _fft->size());
                std::transform(omeg.begin(), omeg.end(), omeg.begin(), [](const auto& value) { return value * 2 * M_PI; });

                const auto& depth = config.receiver_depth();
                const auto nr = depth.points().size();

                types::vector1d_t<types::real_t> tau(nr);
                types::vector2d_t<types::real_t> impulse(nr, types::vector1d_t<types::real_t>(_fft->size()));

                const auto ir = config.reference_index();
                const auto& ref_tau = ir != -1 
                    ? tau[ir] = std::hypot(std::get<0>(depth.points(ir)), std::get<1>(depth.points(ir))) / cm 
                    : 0;

                for (size_t i = 0; i < _impulse_result->size(); ++i) {
                    if (ir != i) {
                        const auto d = (tau[i] = std::hypot(std::get<0>(depth.points(i)), std::get<1>(depth.points(i))) / cm) - ref_tau;

                        auto zip = feniks::zip((*_impulse_result)[i], omeg);
                        std::transform(zip.begin(), zip.end(), _fft->backward_data(),
                            [&](const auto& v) {
                                return std::conj(std::get<0>(v) * std::exp(-1i * d * std::get<1>(v)));
                            }
                        );
                    } else
                        std::transform((*_impulse_result)[i].begin(), (*_impulse_result)[i].end(), _fft->backward_data(), 
                            [](const auto& v) { return std::conj(v); });

                    _fft->execute_backward().normalize_forward();
                    std::memcpy(impulse[i].data(), _fft->forward_data(), _fft->size() * sizeof(types::real_t));
                }

                _owner._meta["tau"] = tau;

                save_impulse(_owner._add_extension(_owner.output / "impulse"), impulse, _owner.binary);

                delete _ix;
                delete _iy;
                delete _impulse_result;
                delete _fft;
            }
        }

    private:

        types::real_t _max;
        jobs_config& _owner;
        types::complex_t _s;
        types::vector1d_t<size_t>* _ix = nullptr;
        types::vector1d_t<types::real_t>* _iy = nullptr;
        types::vector2d_t<types::real_t>* _sel_result = nullptr;
        types::vector2d_t<types::complex_t>* _impulse_result = nullptr;
        acstc::utils::real_fft<types::real_t, types::complex_t>* _fft = nullptr;

        template<typename K0, typename P0, typename KJ>
        auto _perform_init(const K0& k0, const P0& phi_s, const KJ& k_j) {
            const auto init = get_initial_conditions(k0, phi_s, k_j);

            if (_owner.jobs.has_job("init"))
                write_conditions(init, W(_owner._get_filename("init")));

            return init;
        }

        template<typename K0, typename P0>
        auto _perform_modes(const K0& k0, const P0& phi_s) {
            const size_t nm = k0.size();
            _owner._meta["n_modes"].push_back(nm);

            const auto start = std::chrono::system_clock::now();
            auto [k_j, phi_j] = M::make(nm, verbose(2));
            const auto end = std::chrono::system_clock::now();
            verboseln_lv(1, "Modes computing time: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");

            if (k_j.size() > nm) {
                k_j.erase_last(k_j.size() - nm);
                phi_j.erase_last(phi_j.size() - nm);
            }

            if (_owner.jobs.has_job("modes"))
                write_modes(k_j, phi_j, 
                    [writer=W(_owner._get_filename("modes"))](const auto& data) mutable {
                        writer.write(
                            reinterpret_cast<const types::real_t*>(data.data()), 
                            data.size() * sizeof(data[0]) / sizeof(types::real_t)
                        );
                    }
                );

            _owner._meta["k0"].push_back(k0);
            _owner._meta["phi_s"].push_back(phi_s);

            return std::make_tuple(std::move(k_j), std::move(phi_j));
        }        

        template<typename KJ, typename PJ>
        void _perform_rays(const KJ& k_j, const PJ& phi_j) {
            if (_owner.jobs.has_job("rays"))
                save_rays(_owner._get_filename("rays"), _owner.binary, k_j, _owner.step);
        }

        template<typename I, typename K0, typename KJ, typename PJ>
        void _perform_sel(const I& init, const K0& k0, const KJ& k_j, const PJ& phi_j) {
            if (_owner.jobs.has_job("sel"))
                _perform_impulse(init, k0, k_j, phi_j,
                    acstc::utils::ekc_callback(_owner.step,
                        [&sel_result=*_sel_result, i=0](const auto& x, const auto& data) mutable {
                            for (auto [s, v] : feniks::zip(sel_result[i], data))
                                s += std::pow(std::abs(v), 2);
                        }
                    )
                );
            else
                _perform_impulse(init, k0, k_j, phi_j, acstc::utils::nothing_callback());
        }

        template<typename I, typename K0, typename KJ, typename PJ, typename... C>
        void _perform_impulse(const I& init, const K0& k0, const KJ& k_j, const PJ& phi_j, C&&... callbacks) {
            if (_owner.jobs.has_job("impulse"))
                _perform_solution(init, k0, k_j, phi_j, std::forward<C>(callbacks)...,
                    [
                        &depth=config.receiver_depth(), 
                        &result=*_impulse_result, 
                        &ix=*_ix,
                        &iy=*_iy,
                        &fi=config.f_index(),
                        &s=(_s = _fft->backward_data(config.f_index())),
                        last=types::vector1d_t<types::complex_t>(),
                        last_x=config.x0(),
                        li=0,
                        ir=config.reference_index()
                    ](const auto& mx, const auto& data) mutable {
                        const auto x = mx + config.x0();
                        for (; last.size() && li < ix.size() && std::get<0>(depth.points(ix[li])) <= x; ++li) {
                            const auto& [px, py] = depth.points(ix[li]);
                            const auto& [ya, yb] = acstc::utils::find_indices(iy, py);

                            auto& r = result[ix[li]][fi] = acstc::utils::__impl::linear_interpolation::field_point(
                                    last[ya], last[yb], data[ya], data[yb], last_x, x, iy[ya], iy[yb], px, py);

                            if (ix[li] == ir) {
                                const auto buff = s;
                                s = std::abs(r) ? s / r : types::complex_t(0);
                                r = buff;
                            }
                        }
                        last_x = x;
                        last.assign(data.begin(), data.end());
                    }
                );
            else 
                _perform_solution(init, k0, k_j, phi_j, std::forward<C>(callbacks)...);
        }

        template<typename I, typename K0, typename KJ, typename PJ, typename... C>
        void _perform_solution(const I& init, const K0& k0, const KJ& k_j, const PJ& phi_j, C&&... callbacks) {
            if (_owner.jobs.has_job("solution")) {
                W writer(_owner._get_filename("solution"));
                _perform_solve(init, k0, k_j, phi_j, std::forward<C>(callbacks)...,
                    acstc::utils::ekc_callback(_owner.step,
                        [&writer](const auto& x, const auto& data) mutable {
                            writer.write(reinterpret_cast<const types::real_t*>(data.data()), data.size() * 2);
                        }
                    )
                );
            }
            else
                _perform_solve(init, k0, k_j, phi_j, std::forward<C>(callbacks)...);
        }

        template<typename I, typename K0, typename KJ, typename PJ, typename... C>
        void _perform_solve(const I& init, const K0& k0, const KJ& k_j, const PJ& phi_j, C&&... callbacks) {
            if (!_owner.jobs.has_job("solution") && !_owner.jobs.has_job("impulse"))
                return;

            acstc::solver solver(config);

            auto callback = acstc::utils::callbacks(
                acstc::utils::progress_bar_callback(config.nx(), "Solution", verbose(2)),
                std::forward<C>(callbacks)...
            );

            const auto start = std::chrono::system_clock::now();
            solve(solver, init, k0, k_j, phi_j, callback, _owner.num_workers, _owner.buff_size);
            const auto end = std::chrono::system_clock::now();
            verboseln_lv(1, "Elapsed time: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");
        }

    };

};

void validate(boost::any& v, 
              const std::vector<std::string>& values,
              jobs*, int) {
    if (!values.size())
        throw po::validation_error(po::validation_error::at_least_one_value_required);

    if (v.empty()) {
        ::jobs jobs;
        for (const auto& it : values)
            jobs.add_job(it);
        v = boost::any(jobs);
    } else {
        ::jobs& jobs = boost::any_cast<::jobs&>(v);
        for (const auto& it : values)
            jobs.add_job(it);
    }
}

int main(int argc, char* argv[]) {
    try {
        po::positional_options_description positional;
        positional.add("jobs", -1);

        jobs_config jobs_config;
        po::options_description generic("Generic options");
        generic.add_options()
            ("help,h", "Print this message")
            ("verbosity,v", po::value(&acstc::utils::verbosity::instance().level)->default_value(0), "Verbosity level")
            ("config,c", po::value(&jobs_config.config_path)->default_value("config.json"), "Config filename");

        po::options_description output("Output options");
        size_t step;
        std::string output_filename;
        output.add_options()
            ("output,o", po::value(&jobs_config.output)->default_value("output"), "Output filename")
            ("step,s", po::value(&jobs_config.step)->default_value(100)->value_name("k"), "Output every k-th computed row")
            ("binary", "Use binary output");

        po::options_description computation("Computation options");
        size_t num_workers, buff_size;
        computation.add_options()
            ("workers,w", po::value(&jobs_config.num_workers)->default_value(1), "Number of workers for computation")
            ("buff,b", po::value(&jobs_config.buff_size)->default_value(100), "Buff size to be used during multithreaded computation");

        po::options_description options;
        options.add_options()
            ("jobs", po::value(&jobs_config.jobs)->multitoken()->default_value({ "solution" }, "solution"));
        options.add(generic).add(output).add(computation);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).positional(positional).options(options).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            po::options_description desc;
            desc.add(generic).add(output).add(computation);
            std::cout << "Usage: [" << acstc::utils::join(available_jobs.begin(), available_jobs.end(), "|") << 
                "] (=solution) [options]\n" << desc << std::endl;
            return 0;
        }

        jobs_config.binary = vm.count("binary");
        jobs_config.perform();

        return 0;
    }
    catch (const std::exception& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
}
