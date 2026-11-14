#include "cli/orbital_output.hpp"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace kumo_cli {
using chem::kumo::Molecule;
using chem::kumo::Spin;

namespace {
std::size_t number(const std::string& s) {
    if (s.empty() ||
        !std::all_of(s.begin(), s.end(), [](unsigned char c) { return std::isdigit(c); }))
        throw std::invalid_argument("invalid orbital index: " + s);
    std::size_t n = 0;
    const auto value = std::stoull(s, &n);
    if (n != s.size() || value == 0 || value > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("invalid orbital index: " + s);
    return static_cast<std::size_t>(value);
}

// Resolve an absolute MO index or a HOMO/LUMO-relative index (1-based).
std::size_t endpoint(const std::string& s, const Molecule& mol, Spin channel) {
    if (s.empty()) throw std::invalid_argument("empty orbital selector");
    if (std::isdigit(static_cast<unsigned char>(s[0]))) return number(s);
    std::size_t base = 0, offset = 0;
    if (s.rfind("homo", 0) == 0) {
        base = static_cast<std::size_t>(mol.occupied_count(channel));
        offset = 4;
    } else if (s.rfind("lumo", 0) == 0) {
        base = static_cast<std::size_t>(mol.occupied_count(channel)) + 1;
        offset = 4;
    } else
        throw std::invalid_argument("unknown MO selector: " + s);
    if (offset == s.size()) return base;
    const char sign = s[offset];
    if (sign != '+' && sign != '-') throw std::invalid_argument("invalid MO offset: " + s);
    const auto step = number(s.substr(offset + 1));
    if (sign == '-' && step >= base) throw std::out_of_range("MO offset below first orbital: " + s);
    if (sign == '+') {
        if (step > std::numeric_limits<std::size_t>::max() - base)
            throw std::out_of_range("MO index overflow");
        return base + step;
    }
    return base - step;
}

// Expand inclusive ranges and de-duplicate MO selections in user order.
std::vector<std::size_t> select_mos(const std::string& selector, const Molecule& mol,
                                    Spin channel) {
    std::vector<std::size_t> result;
    std::set<std::size_t> seen;
    std::string expanded = selector == "frontier" ? "homo,lumo" : selector;
    std::size_t start = 0;
    while (start < expanded.size()) {
        const auto comma = expanded.find(',', start);
        const auto token =
            expanded.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        const auto colon = token.find(':');
        auto first = endpoint(token.substr(0, colon), mol, channel);
        auto last =
            colon == std::string::npos ? first : endpoint(token.substr(colon + 1), mol, channel);
        if (first > last || last > mol.orbital_count(channel))
            throw std::out_of_range("MO selector outside available orbitals: " + token);
        if (last - first > 10000) throw std::invalid_argument("MO range too large");
        for (auto n = first; n <= last; ++n)
            if (seen.insert(n).second) result.push_back(n);
        if (comma == std::string::npos) break;
        start = comma + 1;
        if (start == expanded.size()) throw std::invalid_argument("trailing comma in MO selector");
    }
    if (result.empty()) throw std::invalid_argument("empty MO selector");
    return result;
}

std::string spin_name(Spin s) {
    return s == Spin::Alpha ? "alpha" : "beta";
}

void list_mos(const Molecule& mol, Spin channel) {
    const auto count = mol.orbital_count(channel);
    const auto& energies = mol.orbital_energies(channel);
    const auto occupied = mol.occupied_count(channel);
    const bool restricted = !mol.unrestricted();
    std::cout << "Spin: " << (restricted ? "restricted" : spin_name(channel))
              << "  orbitals: " << count << "  occupied: " << occupied << "\n";
    std::cout << " MO   Occupation   Energy (hartree)   Label\n";
    for (std::size_t i = 1; i <= count; ++i) {
        std::string label;
        if (i == static_cast<std::size_t>(occupied)) label = "HOMO";
        if (i == static_cast<std::size_t>(occupied) + 1) label = "LUMO";
        std::cout << std::setw(4) << i << "   "
                  << (i <= static_cast<std::size_t>(occupied) ? (restricted ? 2 : 1) : 0)
                  << "            ";
        if (i <= energies.size())
            std::cout << std::setw(15) << std::fixed << std::setprecision(8) << energies[i - 1];
        else
            std::cout << std::setw(15) << "n/a";
        std::cout << "   " << label << "\n";
    }
}

void write_orbitals(const Molecule& mol, const chem::kumo::Grid& grid, const Options& o,
                    Spin channel, Console& console) {
    const auto indices = select_mos(o.selector, mol, channel);
    std::vector<std::string> numbers;
    for (auto n : indices)
        numbers.push_back(std::to_string(n));
    console.begin("Evaluating " + spin_name(channel) + " molecular orbital" +
                  (indices.size() > 1 ? "s " : " ") + join(numbers));
    auto fields = mol.orbitals(grid, indices, channel, console.progress());
    console.end("Evaluated in");
    console.begin("Writing " + std::to_string(indices.size()) + " Cube file" +
                  (indices.size() > 1 ? "s" : ""));
    for (std::size_t j = 0; j < indices.size(); ++j) {
        const auto n = indices[j];
        const auto filename =
            o.prefix + "_" + spin_name(channel) + "_mo_" + std::to_string(n) + ".cube";
        fields[j].write_orbital_cube(filename, n,
                                     "Molecular orbital amplitude (signed), " + spin_name(channel) +
                                         " MO " + std::to_string(n));
        std::cout << filename << "\n";
    }
    console.end("Written in");
}
} // namespace

void list_orbitals(const Molecule& mol, const Options& o) {
    if (!(o.spin.empty() || o.spin == "alpha" || o.spin == "beta" || o.spin == "both"))
        throw std::invalid_argument("--list-mo accepts alpha, beta, or both");
    if (o.spin.empty() || o.spin == "alpha" || o.spin == "both") list_mos(mol, Spin::Alpha);
    if (o.spin == "beta" || o.spin == "both") list_mos(mol, Spin::Beta);
}

void write_orbital_cubes(const Molecule& mol, const chem::kumo::Grid& grid, const Options& o,
                         Console& console) {
    if (!(o.spin.empty() || o.spin == "alpha" || o.spin == "beta" || o.spin == "both"))
        throw std::invalid_argument("--mo accepts --spin alpha, beta, or both");
    if (o.spin.empty() || o.spin == "alpha" || o.spin == "both")
        write_orbitals(mol, grid, o, Spin::Alpha, console);
    if (o.spin == "beta" || o.spin == "both") write_orbitals(mol, grid, o, Spin::Beta, console);
}
} // namespace kumo_cli
