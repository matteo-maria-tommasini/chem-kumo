#include "cli/sampling_output.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace kumo_cli {
using chem::kumo::Molecule;

namespace {
std::vector<std::string> split(const std::string& text, char separator) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    while (true) {
        const auto end = text.find(separator, start);
        parts.push_back(
            text.substr(start, end == std::string::npos ? std::string::npos : end - start));
        if (end == std::string::npos) return parts;
        start = end + 1;
    }
}

// A location given as X,Y,Z in bohr or as @N, the position of atom N (1-based).
Eigen::Vector3d location(const std::string& spec, const Molecule& mol) {
    const auto& atoms = mol.atoms();
    if (!spec.empty() && spec[0] == '@') {
        const auto digits = spec.substr(1);
        if (digits.empty() || digits.size() > 9 ||
            !std::all_of(digits.begin(), digits.end(),
                         [](unsigned char c) { return std::isdigit(c); }))
            throw std::invalid_argument("invalid atom reference '" + spec + "' (use @N)");
        const auto n = std::stoul(digits);
        if (n == 0 || n > atoms.size())
            throw std::invalid_argument("atom " + spec + " does not exist (atoms are @1 to @" +
                                        std::to_string(atoms.size()) + ")");
        return atoms[n - 1].position;
    }
    const auto parts = split(spec, ',');
    Eigen::Vector3d point;
    try {
        if (parts.size() != 3) throw std::invalid_argument("three coordinates");
        for (int k = 0; k < 3; ++k)
            point[k] = real(parts[static_cast<std::size_t>(k)]);
    } catch (const std::exception&) {
        throw std::invalid_argument("invalid location '" + spec +
                                    "' (use X,Y,Z in bohr or @N for atom N)");
    }
    return point;
}

// The locations of --line A:B or --plane A:B:C.
std::vector<Eigen::Vector3d> locations(const std::string& option, const std::string& spec,
                                       std::size_t count, const Molecule& mol) {
    const auto parts = split(spec, ':');
    if (parts.size() != count)
        throw std::invalid_argument(option + " needs " + std::to_string(count) +
                                    " locations separated by ':', for example " +
                                    (count == 2 ? "@1:@2" : "@1:@2:@3"));
    std::vector<Eigen::Vector3d> points;
    for (const auto& part : parts)
        points.push_back(location(part, mol));
    return points;
}

std::string coordinates(const Eigen::Vector3d& p) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << "(" << p.x() + 0.0 << ", " << p.y() + 0.0 << ", "
        << p.z() + 0.0 << ")";
    return out.str();
}

// A table of density properties: '#' comment lines describing the sampling
// and each column, a '#' line with the column names, then one row per point
// with the coordinates (bohr, fixed) and the values (scientific, 10 digits).
class Table {
public:
    Table(const DensitySelection& selection, std::vector<std::string> coordinate_columns)
        : selection_(selection), coordinate_columns_(std::move(coordinate_columns)) {
        for (const auto& [component, suffix] : selection.components)
            for (const auto& property : selection.properties)
                request_.add(density_property(component, property));
    }

    const chem::kumo::Request& request() const noexcept { return request_; }

    void header(std::ostream& out, const std::vector<std::string>& description,
                chem::kumo::Backend used, double esp_tolerance) const {
        for (const auto& line : description)
            out << "# " << line << "\n";
        std::string names;
        for (const auto& name : coordinate_columns_)
            names += (names.empty() ? "" : " ") + name;
        for (const auto& [component, suffix] : selection_.components)
            for (const auto& property : selection_.properties) {
                const auto name = column_name(component, property);
                out << "# " << name << ": " << property.description << ", "
                    << component_label(component) << ", "
                    << method_label(property, used, esp_tolerance) << "\n";
                names += " " + name;
            }
        out << "# " << names << "\n";
    }

    // The value columns of `results`, one vector per column.
    std::vector<std::vector<double>> columns(const chem::kumo::Results& results) const {
        std::vector<std::vector<double>> values;
        for (const auto& [component, suffix] : selection_.components)
            for (const auto& property : selection_.properties)
                values.push_back(scalar_values(results.get(density_property(component, property)),
                                               property.kind));
        return values;
    }

    // Row `q` of `columns`.
    static std::vector<double> at(const std::vector<std::vector<double>>& columns, std::size_t q) {
        std::vector<double> row;
        for (const auto& column : columns)
            row.push_back(column[q]);
        return row;
    }

    static void row(std::ostream& out, const std::vector<double>& coordinates,
                    const std::vector<double>& values) {
        std::ostringstream line;
        line << std::fixed << std::setprecision(8);
        // + 0.0 prints -0 as 0.
        for (std::size_t c = 0; c < coordinates.size(); ++c)
            line << (c > 0 ? " " : "") << std::setw(13) << coordinates[c] + 0.0;
        line << std::scientific << std::setprecision(10);
        for (const double value : values)
            line << " " << std::setw(17) << value;
        out << line.str() << "\n";
    }

private:
    const DensitySelection& selection_;
    std::vector<std::string> coordinate_columns_;
    chem::kumo::Request request_;
};

std::string sampling_summary(const DensitySelection& selection, std::size_t points) {
    return evaluation_summary(selection.components, selection.properties) + " at " +
           Console::grouped(points) + " point" + (points > 1 ? "s" : "");
}

} // namespace

void print_points(const Molecule& mol, const Options& o, const DensitySelection& selection,
                  Console& console) {
    std::vector<Eigen::Vector3d> points;
    for (const auto& spec : o.points)
        points.push_back(location(spec, mol));
    const Table table(selection, {"x", "y", "z"});
    console.begin(sampling_summary(selection, points.size()));
    Eigen::Matrix<double, 3, Eigen::Dynamic> columns(3, static_cast<Eigen::Index>(points.size()));
    for (std::size_t p = 0; p < points.size(); ++p)
        columns.col(static_cast<Eigen::Index>(p)) = points[p];
    const auto results = chem::kumo::Evaluator(mol, evaluation_options(selection, o, console))
                             .evaluate(columns, table.request());
    console.end("Evaluated in");
    const auto values = table.columns(results);
    table.header(std::cout,
                 {"kumo " CHEM_KUMO_VERSION ", " + o.input.filename().string() + ": " +
                      std::to_string(points.size()) + " point" + (points.size() > 1 ? "s" : ""),
                  "x y z: position (bohr)"},
                 results.backend_used(), o.esp_tolerance);
    for (std::size_t p = 0; p < points.size(); ++p)
        Table::row(std::cout, {points[p].x(), points[p].y(), points[p].z()}, Table::at(values, p));
}

// Line rows give t, the signed distance from A along A->B; plane rows give u
// and v, the coordinates in the plane from A along the grid axes (u along
// A->B), with a blank line after each row of constant u, as gnuplot expects
// for surfaces.
void write_profile(const Molecule& mol, const Options& o, const DensitySelection& selection,
                   Console& console) {
    const bool is_line = !o.line.empty();
    const auto ends =
        is_line ? locations("--line", o.line, 2, mol) : locations("--plane", o.plane, 3, mol);
    const auto grid = is_line
                          ? chem::kumo::Grid::line(ends[0], ends[1], o.margin, o.spacing)
                          : chem::kumo::Grid::plane(ends[0], ends[1], ends[2], o.margin, o.spacing);
    const auto nu = static_cast<std::size_t>(grid.counts[0]);
    const auto nv = static_cast<std::size_t>(grid.counts[1]);
    const Eigen::Vector3d u = grid.steps.col(0).normalized(), v = grid.steps.col(1).normalized();

    std::vector<std::string> description;
    std::ostringstream geometry;
    geometry << std::setprecision(6);
    if (is_line)
        geometry << "Line from A " << coordinates(ends[0]) << " to B " << coordinates(ends[1])
                 << " bohr, extended by " << o.margin << " bohr: " << nu << " points "
                 << grid.steps.col(0).norm() << " bohr apart";
    else
        geometry << "Plane of A " << coordinates(ends[0]) << ", B " << coordinates(ends[1])
                 << " and C " << coordinates(ends[2]) << " bohr, margin " << o.margin
                 << " bohr: " << nu << " x " << nv << " points " << grid.steps.col(0).norm()
                 << " x " << grid.steps.col(1).norm() << " bohr apart";
    console.say(geometry.str());
    description.push_back("kumo " CHEM_KUMO_VERSION ", " + o.input.filename().string());
    description.push_back(geometry.str());
    if (is_line)
        description.push_back("t: distance from A towards B (bohr); x y z: position (bohr)");
    else {
        std::ostringstream axes;
        axes << std::setprecision(6) << "u, v: in-plane coordinates from A (bohr), u along A->B; "
             << "B at (" << (ends[1] - ends[0]).dot(u) << ", 0), C at ("
             << (ends[2] - ends[0]).dot(u) << ", " << (ends[2] - ends[0]).dot(v)
             << "); unit normal " << coordinates(grid.steps.col(2)) << "; x y z: position (bohr)";
        description.push_back(axes.str());
    }

    const Table table(selection, is_line ? std::vector<std::string>{"t", "x", "y", "z"}
                                         : std::vector<std::string>{"u", "v", "x", "y", "z"});
    console.begin(sampling_summary(selection, nu * nv));
    const auto results = chem::kumo::Evaluator(mol, evaluation_options(selection, o, console))
                             .evaluate(grid, table.request());
    console.end("Evaluated in");

    const auto filename = o.prefix + (is_line ? "_line.dat" : "_plane.dat");
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("cannot write " + filename);
    table.header(out, description, results.backend_used(), o.esp_tolerance);
    const auto values = table.columns(results);
    for (std::size_t i = 0; i < nu; ++i) {
        if (i > 0 && !is_line) out << "\n";
        for (std::size_t j = 0; j < nv; ++j) {
            const Eigen::Vector3d r = grid.origin + grid.steps.col(0) * static_cast<double>(i) +
                                      grid.steps.col(1) * static_cast<double>(j);
            const Eigen::Vector3d d = r - ends[0];
            std::vector<double> coordinates_row;
            if (is_line)
                coordinates_row = {d.dot(u), r.x(), r.y(), r.z()};
            else
                coordinates_row = {d.dot(u), d.dot(v), r.x(), r.y(), r.z()};
            Table::row(out, coordinates_row, Table::at(values, i * nv + j));
        }
    }
    out.flush();
    if (!out) throw std::runtime_error("write failed: " + filename);
    std::cout << filename << "\n";
}
} // namespace kumo_cli
