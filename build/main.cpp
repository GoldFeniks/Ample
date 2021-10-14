#ifdef _WIN32
#define NOMINMAX
#endif

#include <set>
#include <cmath>
#include <tuple>
#include <chrono>
#include <string>
#include <complex>
#include <cstring>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <iomanip>
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
#include "boundary_conditions.hpp"
#include "utils/multi_optional.hpp"
#include "boost/program_options.hpp"

namespace types = ample::types;
namespace po = boost::program_options;

using nlohmann::json;
using ample::utils::verboseln;
using ample::utils::verboseln_lv;
using namespace std::complex_literals;

static constexpr size_t max_elements = 5;

ample::config<types::real_t> config;

enum class field_group : uint64_t {
    nothing = 0,
    modes = 1 << 1,
    solver = 1 << 2,
    rays = 1 << 3,
    initial = 1 << 4
};

auto operator|(const field_group& a, const field_group& b) {
    return static_cast<field_group>(static_cast<uint64_t>(a) | static_cast<uint64_t>(b));
}

bool operator&(const field_group& a, const field_group& b) {
    return static_cast<uint64_t>(a) & static_cast<uint64_t>(b);
}

inline bool verbose(const size_t& level) {
    return ample::utils::verbosity::instance().level >= level;
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

template<typename T>
struct strings {

    static auto get(const T& values) {
        types::vector1d_t<std::string> result;

        for (const auto& it : values)
            result.push_back(helper.to_string(it));

        return result;
    }

};

template<typename... T>
struct strings<std::tuple<T...>> {

    static auto get(const std::tuple<T...>& values) {
        return get(values, std::make_index_sequence<sizeof...(T)>{});
    }

    template<size_t... I>
    static auto get(const std::tuple<T...>& values, std::integer_sequence<size_t, I...>) {
        types::vector1d_t<std::string> result(sizeof...(I));
        ((result[I] = helper.to_string(std::get<I>(values))), ...);
        return result;
    }

};

template<typename T>
struct strings<types::point<T>> {

    static auto get(const types::point<T>& point) {
        return types::vector1d_t<std::string>{
            helper.to_string(point.x),
            helper.to_string(point.y),
            helper.to_string(point.z)
        };
    }
};

template<typename V>
auto get_strings(const V& values) {
    return strings<V>::get(values);
}

template<typename X, typename Y, typename V>
void print_table(const X& x, const Y& y, const V& values, std::stringstream& stream, const bool& is_torn = false) {
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

    const auto length = std::accumulate(widths.begin(), widths.end(), static_cast<size_t>(0)) + 3 * widths.size() + 1;

    auto* const sep_line = new char[length + 1];
    sep_line[length] = 0;

    std::memset(sep_line, '-', length);

    auto* current = sep_line;
    for (const auto& it : widths) {
        *current = '+';
        current += static_cast<ptrdiff_t>(it + 3);
    }
    *current = '+';

    stream << "        " << sep_line << '\n';

    auto* const buffer = new char[length + 1];
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

    if (is_torn) {
        std::memset(sep_line, '~', length);

        current = sep_line;
        for (const auto& it : widths) {
            *current = '+';
            current += static_cast<ptrdiff_t>(it + 3);
        }
        *current = '+';
    }

    stream << "        " << sep_line << '\n';

    delete[] buffer;
    delete[] sep_line;
}

template<typename T>
void print_values(const char* name, const T& values, std::stringstream& stream) {
    stream << "    " << name << ": [ " << helper.to_string(values[0]);

    for (size_t i = 1; i < std::min(values.size(), max_elements); ++i)
        stream << ", " << helper.to_string(values[i]);

    if (values.size() > max_elements)
        stream << ", ...";

    stream << " ];\n";
}

#define PRINT_FIELD(title, field, stream) \
    if (config.has_##field())             \
        stream << "    " title ": " << helper.to_string(config.field()) << ";\n";

#define PRINT_VALUES(title, field, stream) \
    if (config.has_##field())              \
        print_values(title, config.field(), stream);

#define PRINT_MESH_SPEC(title, field, stream)                                                \
    stream << "    " << title << ": " << helper.to_string(config.field ## 0()) << " <-- "    \
           << (config.has_n##field()                                                          \
              ? helper.to_string(config.n##field())                                          \
              : "")                                                                          \
           << " --> " << helper.to_string(config.field ## 1()) << ";\n";


void print_modes(std::stringstream& stream) {
    stream << "Modal parameters.\n";

    PRINT_FIELD("Mode subset", mode_subset, stream);
    PRINT_FIELD("Points per meter over z", ppm, stream);
    PRINT_FIELD("Richardson extrapolation order", ord_rich, stream);

    PRINT_VALUES("Frequencies, Hz", frequencies, stream);
    PRINT_VALUES("Times, s", times, stream);
    PRINT_VALUES("Source function", source_function, stream);
    PRINT_VALUES("Source spectrum", source_spectrum, stream);

    PRINT_FIELD("Points per meter over z", ppm, stream);
    PRINT_FIELD("Source depth, m", z_s, stream);
    PRINT_FIELD("Complex-valued modes", complex_modes, stream);
    PRINT_FIELD("x-independent modes", const_modes, stream);
    PRINT_FIELD("Use additive depth", additive_depth, stream);

    if (config.has_bottom_c1s() && config.has_bottom_c2s() && config.has_bottom_layers() && config.has_bottom_rhos()) {
        stream << "    Bottom layers (top_speed -> bottom_speed; depth; density):\n";
        for (const auto [c1, c2, z, r] : feniks::zip(
            config.bottom_c1s(),
            config.bottom_c2s(),
            config.bottom_layers(),
            config.bottom_rhos()
        ))
            stream << "        " << helper.to_string(c1) << " -> " << 
                                    helper.to_string(c2) << "; " << 
                                    helper.to_string( z) << "; " << 
                                    helper.to_string( r) << ";\n";
    }

    PRINT_FIELD("Number of water layers", n_layers, stream);
    PRINT_VALUES("Beta parameters", betas, stream);

    const auto nm = config.max_mode();
    stream << "    Maximal number of modes: " << (nm == static_cast<size_t>(-1) ? "All" : helper.to_string(nm)) << ";\n";

    const auto mn = config.n_modes();
    stream << "    Required number of modes: " << (mn == static_cast<size_t>(-1) ? "All" : helper.to_string(mn)) << ";\n";

    if (config.has_bathymetry())
        if (config.const_modes())
            stream << "    Number of points over y: " << config.mny() << ";\n";
        else
            stream << "    Number of points over x: " << config.mnx() << ";\n    Number of points over y: " << config.mny() << ";\n";

    stream << '\n';
}

void print_tapering(std::stringstream& stream) {
    const auto& description = config.tapering();
    const auto vl = THIS_OR_THAT(description.parameters(), types::real_t, "left", "value");
    const auto vr = THIS_OR_THAT(description.parameters(), types::real_t, "right", "value");
    stream << "    Tapering:\n        Type: " << description.type() << ";\n" <<
        "        Left: "  << vl << ";\n" <<
        "        Right: " << vr << ";\n";
}

void print_initial_conditions(std::stringstream& stream) {
    stream << "Initial conditions parameters.\n    Type: ";
    const auto type = config.init();

    if (type == "greene")
        stream << "Greene source";
    else if (type == "gauss")
        stream << "Gaussian source";
    else if (type == "ray_simple") {
        stream << "Ray-based source assuming homogeneous medium;\n";
        PRINT_MESH_SPEC("Angle mesh", a, stream);
        print_tapering(stream);
        stream << '\n';
        return;
    }
    else if (type == "ray") {
        stream << "Ray-based source;\n";
        PRINT_MESH_SPEC("Angle mesh", a, stream);
        PRINT_MESH_SPEC("Natural parameter mesh", l, stream);
        print_tapering(stream);
        stream << '\n';
        return;
    }
    else
        throw std::logic_error(std::string("Unknown initial conditions type: ") + type);

    stream << ";\n\n";
}

void print_solver(std::stringstream& stream) {
    stream << "Solver parameters:\n";

    if (config.has_receivers()) {
        const auto& receivers = config.receivers();
        const auto rn = std::min(max_elements, receivers.size());
        stream << "    Receivers:\n";
        print_table(ample::utils::mesh_1d<size_t>(0, rn - 1, rn), types::vector1d_t<char>{'x', 'y', 'z'},
                    receivers, stream, rn < receivers.size()
        );
    }

    PRINT_FIELD("Source z coordinate", y_s, stream);
    PRINT_FIELD("Width of smooth border over edges", border_width, stream);

    PRINT_MESH_SPEC("x mesh", x, stream);
    PRINT_MESH_SPEC("y mesh", y, stream);

    const auto& description = config.coefficients();
    stream << "    Root approximation coefficients:\n" << 
        "        type: " << description.type() << ";\n" <<
        "        n: " << helper.to_string(description.parameters()["n"].get<size_t>()) << ";\n" <<
        "        m: " << helper.to_string(THIS_OR_THAT(description.parameters(), size_t, "n", "m")) << ";\n";

    stream << '\n';
}

void verbose_config_field_group_parameters(const field_group& group) {
    if (!verbose(3))
        return;

    std::stringstream stream;

    if (group & field_group::modes)
        print_modes(stream);

    if (group & field_group::solver)
        print_solver(stream);

    if (group & field_group::initial)
        print_initial_conditions(stream);

    if (group & field_group::rays) {
        stream << "Rays parameters:\n";
        PRINT_MESH_SPEC("Angle mesh", a, stream);
        PRINT_MESH_SPEC("Natural parameter mesh", l, stream);
    }

    std::cout << stream.str();
}

template<typename F>
auto pass_tapering(const F& func) {
    const auto& description = config.tapering();
    const auto& type = description.type();
    const auto  args = description.has_args("left", "right") ? std::make_tuple("left", "right") : std::make_tuple("value", "value");

    if (type == "percentage")
        return func(description.construct<ample::percentage_tapering<types::real_t>, types::real_t, types::real_t>(std::get<0>(args), std::get<1>(args)));

    if (type == "angled")
        return func(description.construct<ample::angled_tapering<types::real_t>, types::real_t, types::real_t>(std::get<0>(args), std::get<1>(args)));

    throw std::runtime_error(std::string("Unknown tapering type: ") + type);
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const KS& k0, const PS& phi_s,
    const ample::utils::linear_interpolated_data_1d<types::real_t>& k_j) {
    return pass_tapering(
        [&](const auto& tapering) {
            return ample::ray_source(config.x0(), 0., config.y_s(), config.l1(), config.nl(),
                                     config.a0(), config.a1(), config.na(), k0, phi_s, k_j, tapering);
        }
    );
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const size_t& nw, const KS& k0, const PS& phi_s) {
    auto [k_j, phi_j] = config.create_const_modes<types::real_t>({ config.z_s() }, nw, config.n_modes(), verbose(2));

    if (k_j.size() > k0.size())
        k_j.erase_last(k_j.size() - k0.size());

    return get_ray_initial_conditions(k0, phi_s, k_j);
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const KS& k0, const PS& phi_s,
    const ample::utils::linear_interpolated_data_1d<types::real_t, types::complex_t>& k_j) {
    const auto& ys = k_j.get<0>();
    types::vector2d_t<types::real_t> new_k_j(k_j.size(), types::vector1d_t<types::real_t>(ys.size()));

    for (size_t j = 0; j < k_j.size(); ++j)
        std::transform(k_j[j].data().begin(), k_j[j].data().end(), new_k_j[j].begin(), [](const auto& v) { return v.real(); });

    return get_ray_initial_conditions(k0, phi_s, ample::utils::linear_interpolated_data_1d<types::real_t>(ys, new_k_j));
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const KS& k0, const PS& phi_s,
    const ample::utils::linear_interpolated_data_2d<types::real_t>& k_j) {
    types::vector2d_t<types::real_t> new_k_j;
    new_k_j.reserve(k_j.size());

    for (size_t j = 0; j < k_j.size(); ++j)
        new_k_j.emplace_back(k_j[j][0].begin(), k_j[j][0].end());

    return get_ray_initial_conditions(k0, phi_s, ample::utils::linear_interpolated_data_1d<types::real_t>(k_j.get<1>(), new_k_j));
}

template<typename KS, typename PS>
auto get_ray_initial_conditions(const KS& k0, const PS& phi_s,
    const ample::utils::linear_interpolated_data_2d<types::real_t, types::complex_t>& k_j) {
    const auto& ys = k_j.get<1>();
    types::vector2d_t<types::real_t> new_k_j(k_j.size(), types::vector1d_t<types::real_t>(ys.size()));

    for (size_t j = 0; j < k_j.size(); ++j)
        std::transform(k_j[j][0].begin(), k_j[j][0].end(), new_k_j[j].begin(), [](const auto& v) { return v.real(); });

    return get_ray_initial_conditions(k0, phi_s, ample::utils::linear_interpolated_data_1d<types::real_t>(ys, new_k_j));
}

template<typename KS, typename PS>
auto get_simple_initial_conditions(const KS& k0, const PS& phi_s) {
    const auto& init = config.init();

    KS ws(k0.size());
    PS as(phi_s.size());
    std::transform(phi_s.begin(), phi_s.end(), as.begin(), [](const auto& phi) { return phi / (2 * std::sqrt(M_PI)); });
    std::transform(k0.begin(), k0.end(), ws.begin(), [](const auto& k0) { return 1. / std::pow(k0, 2); } );

    if (init == "greene")
        return ample::greene_source<types::complex_t>(config.y0(), config.y1(), config.ny(), config.y_s(), as, ws);

    if (init == "gauss")
        return ample::gaussian_source<types::complex_t>(config.y0(), config.y1(), config.ny(), config.y_s(), as, ws);

    if (init == "ray_simple")
        return pass_tapering(
            [&](const auto& tapering) {
                return ample::simple_ray_source(config.x0(), config.a0(), config.a1(), k0, phi_s, tapering);
            }
        );

    throw std::runtime_error(std::string("Unknown initial conditions type: ") + init);
}

template<typename KS, typename PS, typename KJ>
auto get_initial_conditions(const KS& k0, const PS& phi_s, const KJ& k_j) {
    const auto& init = config.init();

//    if (init == "ray")
//        return get_ray_initial_conditions(k0, phi_s, k_j);

    return get_simple_initial_conditions(k0, phi_s);
}

template<typename KS, typename PS>
auto get_initial_conditions(const size_t& nw, const KS& k0, const PS& phi_s) {
    const auto& init = config.init();

//    if (init == "ray")
//        return get_ray_initial_conditions(nw, k0, phi_s);

    return get_simple_initial_conditions(k0, phi_s);
}

template<typename V, typename W>
void write_modes(const ample::utils::linear_interpolated_data_1d<types::real_t, V>& modes, W&& writer) {
    for (size_t i = 0; i < modes.size(); ++i)
        writer(modes[i].data());
}

template<typename V, typename W>
void write_modes(const ample::utils::linear_interpolated_data_2d<types::real_t, V>& modes, W&& writer) {
    for (size_t i = 0; i < modes.size(); ++i)
        for (const auto& it : modes[i].data())
            writer(it);
}

template<typename V, typename W>
void write_modes(const ample::utils::linear_interpolated_data_3d<types::real_t, V>& modes, W&& writer) {
    for (size_t i = 0; i < modes.size(); ++i)
        for (const auto& field : modes[i].data())
            for (const auto& row : field)
                writer(row);
}

template<typename V, typename W>
void write_strided(const V& v, const size_t& k, W&& writer) {
    auto [begin, end] = ample::utils::stride(v.begin(), v.end(), k);
    writer(begin, end);
}

template<typename RX, typename RY, typename W>
void write_rays(const RX& rx, const RY& ry, const size_t& n, const size_t& kr, const size_t& kc, W&& writer) {
    for (const auto [rx, ry] : feniks::zip(rx, ry))
        write_strided(feniks::zip(rx.data(), ry.data()), kr,
            [&writer, &kc](auto begin, const auto& end) mutable {
                while (begin != end) {
                    writer.before_write();
                    const auto& [x, y] = *begin;
                    write_strided(feniks::zip(x, y), kc,
                        [&writer](auto begin, const auto& end) mutable {
                            while (begin != end) {
                                const auto& [a, b] = *begin;
                                writer.write_one(a);
                                writer.write_one(b);
                                ++begin;
                            }
                        }
                    );
                    writer.after_write();
                    ++begin;
                }
            }
        );
}

template<typename V, typename W>
void write_conditions(const V& values, const size_t& k, W&& writer) {
    for (const auto& it : values)
        write_strided(it, k, writer);
}

template<typename W>
void write_impulse(const types::vector2d_t<types::real_t>& impulse, W&& writer) {
    for (const auto& it : impulse)
        writer.write(it);
}

template<typename V, typename W>
void write_vector(const types::vector1d_t<V>& data, W&& writer) {
    writer.write(reinterpret_cast<const types::real_t*>(data.data()), data.size() * sizeof(data[0]) / sizeof(types::real_t));
}

template<typename S, typename K, typename V, typename I, typename C>
void solve(S& solver, const I& init, const K& k0,
           const ample::utils::linear_interpolated_data_1d<types::real_t, V>& k_j,
           const ample::utils::linear_interpolated_data_2d<types::real_t>& phi_j,
           C&& callback, const size_t& num_workers, const size_t& buff_size) {
    solver.solve(init, k0, k_j, phi_j, callback, num_workers, buff_size);
}

template<typename S, typename K, typename V, typename I, typename C>
void solve(S& solver, const I& init, const K& k0,
           const ample::utils::linear_interpolated_data_2d<types::real_t, V>& k_j,
           const ample::utils::linear_interpolated_data_3d<types::real_t>& phi_j,
           C&& callback, const size_t& num_workers, const size_t& buff_size) {
    solver.solve(init, k0, k_j, phi_j, callback, num_workers, buff_size);
}

const std::set<std::string> available_jobs {
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
    jobs(const std::initializer_list<std::string> init) : _jobs(init) {}

    void add_job(const std::string& job) {
        if (available_jobs.find(job) == available_jobs.end())
            throw std::logic_error(std::string("Unknown job type: ") + job);
        _jobs.insert(job);
    }

    [[nodiscard]] bool has_job(const std::string& name) const {
        return _jobs.find(name) != _jobs.end();
    }

    [[nodiscard]] auto size() const {
        return _jobs.size();
    }

    [[nodiscard]] const auto& raw() const {
        return _jobs;
    }

private:

    std::unordered_set<std::string> _jobs;

};

class jobs_config {

public:

    ::jobs jobs;
    bool binary;
    size_t row_step, col_step, num_workers, buff_size;
    std::filesystem::path output, config_path;

    void command_line_arguments(types::vector1d_t<std::string> args) {
        _meta["command_line_arguments"] = std::move(args);
    }

    void perform() {
        config.update_from_file(config_path.generic_string());
        ample::utils::progress_bar::clear_on_end = true;
        std::filesystem::create_directories(output);

        field_group group = field_group::nothing;

        _prep("init", field_group::initial, group);
        _prep("rays", field_group::rays, group);
        _prep("modes", "phi_j", field_group::modes, group);
        _prep("modes", "k_j", field_group::modes, group);
        _prep("modes", "k0", field_group::modes, group);
        _prep("modes", "phi_s", field_group::modes, group);
        _prep("solution", field_group::modes | field_group::solver | field_group::initial, group);

        if (jobs.has_job("impulse"))
            group = group | field_group::modes | field_group::solver | field_group::initial;

        verbose_config_field_group_parameters(group);

        _meta["f"] = json::array();
        _meta["jobs"] = jobs.raw();
        _meta["outputs"] = json::array();
        _meta["original_config_path"] = config_path.generic_string();

        const auto start = std::chrono::system_clock::now();
        _pick_writer();
        const auto end = std::chrono::system_clock::now();

        _meta["computation_time"] = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.;

        const auto dimx = ample::utils::lazy_value<json>([&]() { return _dimension(config.x0(), config.x1(), config.nx(), row_step); });
        const auto dimy = ample::utils::lazy_value<json>([&]() { return _dimension(config.y0(), config.y1(), config.ny(), col_step); });
        const auto dimz = ample::utils::lazy_value<json>([&]() { return _dimension(config.z0(), config.z1(), config.nz()); });
        const auto dimm = ample::utils::lazy_value<json>([&]() { return _dimensions(_n_modes); });

        const auto files = ample::utils::make_vector(_meta["f"].get<types::vector1d_t<types::real_t>>(),
            [this](const auto& value) { return _add_extension(helper.to_string(value)); });

        if (jobs.has_job("sel"))
            _meta["outputs"].push_back(_get_meta_for("sel", { dimx(), dimy(), dimz() }, _add_extension(std::string("sel"))));

        if (jobs.has_job("init"))
            _save_meta_for("init", { dimm(), dimy() }, files);        

        if (jobs.has_job("rays"))
            _save_meta_for("rays", { 
                    dimm(),
                    _dimension(config.a0(), config.a1(), config.na(), row_step),
                    _dimension(config.l0(), config.l1(), config.nl(), col_step)
                }, files
            );

        if (jobs.has_job("modes")) {
            if (config.const_modes()) {
                _save_meta_for("phi_j",
                   {
                       dimm(),
                       _dimension(config.y0(), config.y1(), config.mny()),
                       _dimension(config.z0(), config.z1(), config.mnz())
                   }, files
                );
                _save_meta_for("k_j", config.complex_modes() ? "complex_k_j" : "k_j",
                   {
                       dimm(),
                       _dimension(config.y0(), config.y1(), config.mny())
                   }, files
                );
            } else {
                _save_meta_for("phi_j",
                   {
                       dimm(),
                       _dimension(config.x0(), config.x1(), config.mnx()),
                       _dimension(config.y0(), config.y1(), config.mny()),
                       _dimension(config.z0(), config.z1(), config.mnz())
                   }, files
                );
                _save_meta_for("k_j", config.complex_modes() ? "complex_k_j" : "k_j",
                   {
                       dimm(),
                       _dimension(config.x0(), config.x1(), config.mnx()),
                       _dimension(config.y0(), config.y1(), config.mny())
                   }, files
                );
            }

            _save_meta_for("k0", config.complex_modes() ? "complex_k0" : "k0", { dimm() }, files);
            _save_meta_for("phi_s", { dimm() }, files);
        }

        if (jobs.has_job("impulse"))
            _meta["outputs"].push_back(_get_meta_for("impulse", { 
                    _dimension(config.receivers()),
                    _dimension(config.times().front(), config.times().back(), config.times().size())
                }, _add_extension(std::string("impulse")))
            );

        if (jobs.has_job("solution"))
            _save_meta_for("solution", { _n_modes.size(), dimx(), dimy(), dimz() }, files);

        std::ofstream out(output / "meta.json");
        out << std::setw(4) << _meta << std::endl;

        config.save(output);
    }

private:

    json _meta;
    types::vector1d_t<size_t> _n_modes;

    template<typename V>
    static json _dimensions(const V& values) {
        return ample::utils::make_vector(values, [](const auto& value) { return json{ { "n", value } }; });
    }

    template<typename V>
    static json _dimension(const V& values) {
        return {
            { "n", values.size() },
            { "values", values }
        };
    }

    static json _dimension(const types::real_t& a, const types::real_t& b, const size_t& n) {
        return {
            { "n", n },
            {
                "bounds",
                {
                    { "a", a },
                    { "b", b },
                    { "d", n > 1 ? (b - a) / (n - 1) : size_t(0) }
                }
            }
        };
    }

    static json _dimension(const types::real_t& a, const types::real_t& b, const size_t& n, const size_t& step) {
        const auto h = (b - a) / (n - 1);
        const auto m = (n - 1) / step + 1;

        return _dimension(a, a + h * step * (m - 1), m);
    }

    void _save_meta_for(const std::string& type, const json& dimensions, const types::vector1d_t<std::string>& files) {
        _save_meta_for(type, type, dimensions, files);
    }

    void _save_meta_for(const std::string& path, const std::string& type, const json& dimensions, const types::vector1d_t<std::string>& files) {
        const auto filename = output / path / "meta.json";
        std::ofstream out(filename);

        out << std::setw(4) << _get_meta_for(type, dimensions, files);        

        _meta["outputs"].push_back((filename.parent_path().filename() / filename.filename()).generic_string());
    }

    json _get_meta_for(const std::string& type, const json& dimensions, const std::string& file) {
        return {
            { "type", type },
            { "dimensions", dimensions },
            { "values", file },
            { "binary", binary }
        };
    }

    json _get_meta_for(const std::string& type, const json& dimensions, const types::vector1d_t<std::string>& files) {        
        return {
            { "type", type },
            { "dimensions", dimensions },
            { "values", files },
            { "binary", binary }
        };
    }

    void _make_output_folder(const char* folder) const {
        std::filesystem::create_directories(output / folder);
    }

    void _prep(const char* name, const char* folder, const field_group& params, field_group& group) const {
        if (jobs.has_job(name)) {
            _make_output_folder(folder);
            group = group | params;
        }
    }

    void _prep(const char* name, const field_group& params, field_group& group) const {
        _prep(name, name, params, group);
    }

    template<typename T>
    [[nodiscard]] T _add_extension(T path) const {
        path += binary ? ".bin" : ".txt";
        return path;
    }

    auto _get_filename(const char* name) const {
        return _add_extension(output / name / helper.to_string(config.f()));
    }

    void _pick_writer() {
        if (binary)
            _pick_const<ample::utils::binary_writer>();
        else
            _pick_const<ample::utils::text_writer>();
    }

    template<template<typename> typename W>
    void _pick_const() {
        if (config.const_modes())
            _pick_complex<W, true>();
        else
            _pick_complex<W, false>();
    }

    template<template<typename> typename W, bool Const>
    void _pick_complex() {
        if (config.complex_modes())
            performer<W, make_modes<Const, types::complex_t>>(*this).perform();
        else
            performer<W, make_modes<Const, types::real_t>>(*this).perform();

    }

    template<bool Const, typename T>
    struct make_modes {

        static auto make(const size_t& nw, const size_t& nm, const bool& show_progress) {
            if constexpr (Const)
                return config.create_const_modes<T>(nw, nm, show_progress);
            else
                return config.create_modes<T>(nw, nm, show_progress);
        }

        static auto make_source() {
            return config.create_source_modes<T>(config.n_modes());
        }

    };

    template<template<typename> typename W, typename M>
    class performer {

    public:

        performer(jobs_config& owner) : _owner(owner) {}

        void perform() {
            const auto has_sel = _owner.jobs.has_job("sel");
            const auto has_impulse = _owner.jobs.has_job("impulse");

            if (has_impulse) {
                _load_source_spectrum();

                const auto& receivers = config.receivers();
                const auto nr = receivers.size();

                _impulse_result = new types::vector2d_t<types::complex_t>(nr,
                    types::vector1d_t<types::complex_t>(config.frequencies().size(), static_cast<types::real_t>(0))
                );

                _ix = new types::vector1d_t<size_t>(nr);
                std::iota(_ix->begin(), _ix->end(), 0);
                std::sort(_ix->begin(), _ix->end(),
                    [&receivers](const auto& a, const auto& b) {
                        return receivers[a].x < receivers[b].x;
                    }
                );

                _iy = new types::vector1d_t<types::real_t>(ample::utils::mesh_1d(config.y0(), config.y1(), config.ny()));
                _iz = new types::vector1d_t<types::real_t>(ample::utils::mesh_1d(config.z0(), config.z1(), config.nz()));
            }

            if (has_sel) {
                if (!_source_spectrum.has_value())
                    _load_source_spectrum();

                _sel_result = new types::vector3d_t<types::real_t>(
                        (config.nx() - 1) / _owner.row_step + 1,
                    types::vector2d_t<types::real_t>(
                        (config.ny() - 1) / _owner.col_step + 1,
                        types::vector1d_t<types::real_t>(config.nz(), 0)
                    )
                );

                _sel_buffer = new types::vector3d_t<types::complex_t>(
                        (config.nx() - 1) / _owner.row_step + 1,
                    types::vector2d_t<types::complex_t>(
                        (config.ny() - 1) / _owner.col_step + 1,
                        types::vector1d_t<types::complex_t>(config.nz(), 0)
                    )
                );
            }

            const auto [f0, f1] = config.sel_range();
            ample::utils::progress_bar pbar(config.frequencies().size(), "Frequency", verbose(2), ample::utils::progress_bar::on_end::leave);

            for (const auto& fi : pbar) {
                config.index(fi);

                const auto f = config.f();
                if ((
                        _source_spectrum.has_value() && std::abs(_source_spectrum[fi]) < _max * config.tolerance() || 
                        config.sel_strict()
                    ) && !(config.sel_strict() && f0 <= f && f1 >= f))
                    continue;

                _owner._meta["f"].push_back(f);

                if (_source_spectrum.has_value())
                    _s = _source_spectrum[fi];

                const auto [k0, phi_s] = M::make_source();
                const auto [k_j, phi_j] = _perform_modes(k0, phi_s);

                _perform_rays(k_j, phi_j);

                if (!(has_impulse || has_sel || _owner.jobs.has_job("solution") || _owner.jobs.has_job("init")))
                    continue;
                
                const auto init = _perform_init(k0, phi_s, k_j);

                _perform_sel(init, k0, k_j, phi_j);

                if (k0.empty())
                    continue;

                if (has_impulse) {
                    auto& result = *_impulse_result;
                    for (size_t i = 0; i < result.size(); ++i)
                        if (i != config.reference_index())
                            result[i][fi] *= _s;
                }

                if (has_sel) {
                    for (auto [vyz, byz] : feniks::zip(*_sel_result, *_sel_buffer))
                        for (auto [vz, bz] : feniks::zip(vyz, byz))
                            for (auto [v, b] : feniks::zip(vz, bz))
                                v += std::pow(std::abs(b * _s), 2);
                }
            }

            if (has_sel) {
                const auto size = has_impulse ? _fft->size() : _source_spectrum.size();
                const auto coef = 2 * config.dT() / (size * size);

                for (auto& x : *_sel_result)
                    for (auto& y : x)
                        for (auto& z : y)
                            z *= coef;

                W<types::real_t> writer(_owner._add_extension(_owner.output / "sel"));
                for (const auto& y : *_sel_result)
                    for (const auto& z : y)
                        writer.write(z);

                delete _sel_result;
                delete _sel_buffer;
            }

            if (has_impulse) {
                const auto h  = config.bathymetry().point(0, config.y_s());
                const auto cs = config.hydrology().line(0, 0, h, static_cast<size_t>(h));
                const auto cw = config.c_win();

                auto omeg = ample::utils::mesh_1d(0., 1 / config.dt(), _fft->size());
                std::transform(omeg.begin(), omeg.end(), omeg.begin(), [](const auto& value) { return value * 2 * M_PI; });

                const auto& receivers = config.receivers();
                const auto nr = receivers.size();

                types::vector1d_t<types::real_t> tau(nr);
                types::vector2d_t<types::real_t> impulse(nr, types::vector1d_t<types::real_t>(_fft->size()));

                for (size_t i = 0; i < _impulse_result->size(); ++i) {
                    const auto d = (tau[i] = std::hypot(receivers[i].x, receivers[i].y) / cw);

                    auto zip = feniks::zip((*_impulse_result)[i], omeg);
                    std::transform(zip.begin(), zip.end(), _fft->backward_data(),
                        [&](const auto& v) {
                            return std::conj(std::get<0>(v) * std::exp(-1i * d * std::get<1>(v)));
                        }
                    );

                    _fft->execute_backward().normalize_forward();
                    std::memcpy(impulse[i].data(), _fft->forward_data(), _fft->size() * sizeof(types::real_t));
                }

                _owner._meta["tau"] = tau;

                write_impulse(impulse, W<types::real_t>(_owner._add_extension(_owner.output / "impulse").generic_string()));

                delete _ix;
                delete _iy;
                delete _iz;
                delete _impulse_result;
                delete _fft;
            }
        }

    private:

        types::real_t _max;
        jobs_config& _owner;
        types::complex_t _s;
        types::vector1d_t<size_t>* _ix = nullptr;
        types::vector3d_t<types::real_t>* _sel_result = nullptr;
        ample::utils::span<const types::complex_t> _source_spectrum;
        types::vector1d_t<types::real_t>* _iy = nullptr, *_iz = nullptr;
        ample::utils::real_fft<types::real_t, types::complex_t>* _fft = nullptr;
        types::vector2d_t<types::complex_t>* _impulse_result = nullptr;
        types::vector3d_t<types::complex_t>* _sel_buffer = nullptr;

        void _load_source_spectrum() {
            ample::utils::dynamic_assert(!(config.has_source_function() && config.has_source_spectrum()),
                                         "Only one of source function or spectrum can be provided to compute impulse");

            if (config.has_source_function()) {
                _fft = new ample::utils::real_fft<types::real_t, types::complex_t>(static_cast<int>(config.source_function().size()));
                std::memcpy(_fft->forward_data(), config.source_function().data(), _fft->size() * sizeof(types::real_t));
                _fft->execute_forward();

                std::transform(std::as_const(*_fft).backward_data(), _fft->backward_data_end(), _fft->backward_data(), 
                    [](const auto& v) { 
                        return std::conj(v); 
                    }
                );

                config.frequencies(
                    ample::utils::make_vector_i(static_cast<size_t>(_fft->size()) / 2 + 1,
                        [dT=config.dT()](const size_t& i) {
                            return i / dT;
                        }
                    )
                );
            } else if (config.has_source_spectrum()) {
                _fft = new ample::utils::real_fft<types::real_t, types::complex_t>(static_cast<int>(config.source_spectrum().size()));
                std::memcpy(_fft->backward_data(), config.source_spectrum().data(), _fft->size() * sizeof(types::complex_t));
            } else
                throw std::runtime_error("Either source function or spectrum must be provided");

            _source_spectrum.assign(_fft->backward_data(), config.frequencies().size());

            _max = std::abs(*std::max_element(std::as_const(*_fft).backward_data(), _fft->backward_data_end(),
                [](const auto& a, const auto& b) {
                    return std::abs(a) < std::abs(b);
                }
            ));
        }

        template<typename K0, typename P0, typename KJ>
        auto _perform_init(const K0& k0, const P0& phi_s, const KJ& k_j) {
            const auto init = get_initial_conditions(k0, phi_s, k_j);

            if (_owner.jobs.has_job("init"))
                write_conditions(init.make(config.y0(), config.y1(), config.ny(), k0.size()), _owner.col_step, W<types::complex_t>(_owner._get_filename("init")));

            return init;
        }

        template<typename K0, typename P0>
        auto _perform_modes(const K0& k0, const P0& phi_s) {
            const size_t nm = k0.size();
            _owner._n_modes.push_back(nm);

            using modes_t = decltype(M::make(size_t(0), size_t(0), false));
            std::decay_t<std::tuple_element_t<0, modes_t>> k_j;
            std::decay_t<std::tuple_element_t<1, modes_t>> phi_j;

            if (nm > 0) {
                const auto start = std::chrono::system_clock::now();
                std::tie(k_j, phi_j) = M::make(_owner.num_workers, nm, verbose(2));
                const auto end = std::chrono::system_clock::now();
                verboseln_lv(1, "Modes computing time: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");

                if (k_j.size() > nm) {
                    k_j.erase_last(k_j.size() - nm);
                    phi_j.erase_last(phi_j.size() - nm);
                }
            }

            if (_owner.jobs.has_job("modes")) {
                write_modes(k_j,
                    [writer=W<types::real_t>(_owner._get_filename("k_j"))](const auto &data) mutable {
                        write_vector(data, writer);
                });
                write_modes(phi_j, W<types::real_t>(_owner._get_filename("phi_j")));

                write_vector(k0,    W<types::real_t>(_owner._get_filename("k0")));
                write_vector(phi_s, W<types::real_t>(_owner._get_filename("phi_s")));
            }

            return std::make_tuple(std::move(k_j), std::move(phi_j));
        }        

        template<typename KJ, typename PJ>
        void _perform_rays(const KJ& k_j, const PJ& phi_j) {
            if (_owner.jobs.has_job("rays")) {
                const auto na = config.na();
                const auto nl = config.nl();
                const auto nm = k_j.size();

                auto writer = W<types::real_t>(_owner._get_filename("rays"));
                writer.before_write();

                auto write_rays = [
                    &writer,
                    rs=_owner.row_step,
                    cs=_owner.col_step,
                    p=size_t(0),
                    last_k=size_t(-1)
                ](const auto& j, const auto& i, const auto& k, const auto& x, const auto& y, const auto& a, const auto& l) mutable {
                    if (i % rs != 0)
                        return;

                    if (last_k != k) {
                        last_k = k;
                        p = 0;

                        writer.after_write();
                        writer.before_write();
                    }

                    if (p % cs == 0) {
                        writer.write_one(x);
                        writer.write_one(y);
                    }

                    ++p;
                };

                ample::rays::compute(config.x0(), config.y_s(), config.l1(), nl, config.a0(), config.a1(), na, k_j, write_rays, verbose(2));
            }
        }

        template<typename I, typename K0, typename KJ, typename PJ>
        void _perform_sel(const I& init, const K0& k0, const KJ& k_j, const PJ& phi_j) {
            if (_owner.jobs.has_job("sel"))
                _perform_impulse(init, k0, k_j, phi_j,
                    ample::utils::ekc_callback(_owner.row_step,
                        [&sel_buffer=*_sel_buffer, i=size_t(0), this](const auto& x, const auto& data) mutable {
                            size_t j = 0;

                            auto [begin, end] = ample::utils::stride(data.begin(), data.end(), _owner.col_step);
                            while (begin != end) {
                                sel_buffer[i][j++] = *begin;
                                ++begin;
                            }

                            ++i;
                        }
                    )
                );
            else
                _perform_impulse(init, k0, k_j, phi_j, ample::utils::nothing_callback());
        }

        template<typename I, typename K0, typename KJ, typename PJ, typename... C>
        void _perform_impulse(const I& init, const K0& k0, const KJ& k_j, const PJ& phi_j, C&&... callbacks) {
            if (_owner.jobs.has_job("impulse"))
                _perform_solution(init, k0, k_j, phi_j, std::forward<C>(callbacks)...,
                    [
                        &receivers=config.receivers(),
                        &result=*_impulse_result, 
                        &ix=*_ix,
                        &iy=*_iy,
                        &iz=*_iz,
                        &fi=config.index(),
                        &s=_s,
                        last=types::vector2d_t<types::complex_t>(),
                        last_x=config.x0(),
                        li=0,
                        ir=config.reference_index()
                    ](const auto& x, const auto& data) mutable {
                        for (; !last.empty() && li < ix.size() && receivers[ix[li]].x <= x; ++li) {
                            const auto& [px, py, pz] = receivers[ix[li]];
                            const auto& [ya, yb] = ample::utils::find_indices(iy, py);
                            const auto& [za, zb] = ample::utils::find_indices(iz, pz);

                            auto& r = result[ix[li]][fi] = ample::utils::_impl::linear_interpolation::area_point(
                                    last[ya][za], data[ya][za], last[ya][zb], data[ya][zb],
                                    last[yb][za], data[yb][za], last[yb][zb], data[yb][zb],
                                    last_x, x, iy[ya], iy[yb], iz[za], iz[zb], px, py, pz);

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
                W<types::real_t> writer(_owner._get_filename("solution"));
                _perform_solve(init, k0, k_j, phi_j, std::forward<C>(callbacks)...,
                    ample::utils::ekc_callback(_owner.row_step,
                        [&writer, this](const auto& x, const auto& data) mutable {
                            for (size_t i = 0; i < data.size(); i += _owner.col_step)
                                writer.write(reinterpret_cast<const types::real_t*>(data[i].data()), data[i].size() * 2);
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

            const auto descriptor = config.boundary_conditions();

            ample::solver solver(
                descriptor.construct<
                    ample::pml_boundary_conditions<types::real_t>, size_t, ample::pml_function<types::real_t>>("width" ,"function"),
                config
            );

            auto callback = ample::utils::callbacks(
                ample::utils::progress_bar_callback(config.nx(), "Solution", verbose(2)),
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
    if (values.empty())
        throw po::validation_error(po::validation_error::at_least_one_value_required);

    if (v.empty()) {
        ::jobs jobs;
        for (const auto& it : values)
            jobs.add_job(it);
        v = boost::any(jobs);
    } else {
        auto& jobs = boost::any_cast<::jobs&>(v);
        for (const auto& it : values)
            jobs.add_job(it);
    }
}

#ifdef WIN32
#include "boost/locale.hpp"

int wmain(int argc, const wchar_t* argv[]) {
    types::vector1d_t<std::string> args;
    for (int i = 0; i < argc; ++i)
        args.emplace_back(boost::locale::conv::utf_to_utf<char>(std::wstring(argv[i])));
#else
int main(int argc, const char* argv[]) {
    types::vector1d_t<std::string> args;
    for (size_t i = 0; i < argc; ++i)
        args.emplace_back(std::string(argv[i]));
#endif
    try {
        jobs_config jobs_config;
        jobs_config.command_line_arguments(args);
        args.erase(args.begin());

        po::positional_options_description positional;
        positional.add("jobs", -1);

        po::options_description generic("Generic options");
        generic.add_options()
            ("help,h", "Print this message")
            ("verbosity,v", po::value(&ample::utils::verbosity::instance().level)->default_value(0), "Verbosity level")
            ("config,c", po::value(&jobs_config.config_path)->default_value("config.json"), "Config filename");

        po::options_description output("Output options");
        output.add_options()
            ("output,o", po::value(&jobs_config.output)->default_value("output"), "Output filename")
            ("row_step", po::value(&jobs_config.row_step)->default_value(10)->value_name("k"), "Output every k-th computed row")
            ("col_step", po::value(&jobs_config.col_step)->default_value(1)->value_name("k"), "Output every k-th computed column")
            ("binary", "Use binary output");

        po::options_description computation("Computation options");
        size_t num_workers, buff_size;
        computation.add_options()
            ("workers,w", po::value(&jobs_config.num_workers)->default_value(1), "Number of workers for computation")
            ("buff,b", po::value(&jobs_config.buff_size)->default_value(100), "Buff size to be used during multi-threaded computation");

        po::options_description options;
        options.add_options()
            ("jobs", po::value(&jobs_config.jobs)->multitoken()->default_value({ "solution" }, "solution"));
        options.add(generic).add(output).add(computation);

        po::variables_map vm;
        po::store(po::command_line_parser(args).positional(positional).options(options).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            po::options_description desc;
            desc.add(generic).add(output).add(computation);
            std::cout << "Usage: [ [" << ample::utils::join_it(available_jobs.begin(), available_jobs.end(), "|") <<
                "], ... ] (=solution) [options]\n" << desc << std::endl;
            return 0;
        }

        jobs_config.binary = vm.count("binary");
        jobs_config.perform();

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
