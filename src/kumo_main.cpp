#include "chem/kumo.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using chem::kumo::Molecule;
using chem::kumo::Spin;
namespace fs = std::filesystem;
namespace {
// CLI options only: numerical evaluation stays in libkumo, not in main().
struct Options {
    fs::path input;
    std::string selector;
    std::string spin;
    std::string prefix;
    double spacing = 0.20;
    double margin = 5.0;
    bool mo_mode = false;
    bool list = false;
    bool output_set = false;
};
void usage() {
    std::cout << "kumo " CHEM_KUMO_VERSION " - ground-state electron densities and molecular orbitals\n\n"
      << "Usage: kumo FILE.fchk[.gz] [options]\n\n"
      << "  (default)                  Write total electron density\n"
      << "  --spin total|alpha|beta|difference|all   Density component(s)\n"
      << "  --mo N[,M:K,...]           Write 1-based molecular orbitals\n"
      << "  --mo homo|lumo|frontier    Frontier orbital selection\n"
      << "  --mo homo-2:lumo+2         Frontier-relative range\n"
      << "  --spin alpha|beta|both     MO spin channel (default alpha)\n"
      << "  --list-mo                  List available orbitals and energies\n"
      << "  --spacing BOHR             Grid step (default 0.20 bohr)\n"
      << "  --margin BOHR              Bounding-box margin (default 5.0 bohr)\n"
      << "  --output PREFIX            Output filename prefix\n"
      << "  --help                     Show help\n"
      << "  --version                  Show version\n";
}
std::string lower(std::string s) {
    for (auto& c:s) c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}
double real(const std::string& s) {
    std::size_t n=0;
    const auto v=std::stod(s,&n);
    if (n!=s.size() || !std::isfinite(v)) throw std::invalid_argument("invalid real number: "+s);
    return v;
}
Options parse(int argc,char** argv) {
    Options o;
    for (int i=1;i<argc;++i) {
        const std::string a=argv[i];
        auto value=[&]() -> std::string {
            if (i+1>=argc) throw std::invalid_argument("missing argument after "+a);
            return argv[++i];
        };
        if (a=="--mo") { if (o.mo_mode) throw std::invalid_argument("--mo specified twice"); o.selector=lower(value());o.mo_mode=true; }
        else if (a=="--spin") { if (!o.spin.empty()) throw std::invalid_argument("--spin specified twice");o.spin=lower(value()); }
        else if (a=="--spacing") o.spacing=real(value());
        else if (a=="--margin") o.margin=real(value());
        else if (a=="--output") { o.prefix=value();o.output_set=true; }
        else if (a=="--list-mo") o.list=true;
        else if (!a.empty() && a[0]=='-') throw std::invalid_argument("unknown option: "+a);
        else if (o.input.empty()) o.input=a;
        else throw std::invalid_argument("unexpected positional argument: "+a);
    }
    if (o.input.empty()) throw std::invalid_argument("missing input FCHK file");
    if (o.mo_mode && o.list) throw std::invalid_argument("choose --mo or --list-mo");
    if (!(o.spacing>0) || o.margin<0) throw std::invalid_argument("invalid grid spacing or margin");
    if (!o.output_set) o.prefix=o.input.stem().string();
    if (o.prefix.empty()) throw std::invalid_argument("empty output prefix");
    return o;
}
std::size_t number(const std::string& s) {
    if (s.empty() || !std::all_of(s.begin(),s.end(),[](unsigned char c){return std::isdigit(c);}))
        throw std::invalid_argument("invalid orbital index: "+s);
    std::size_t n=0;
    const auto value=std::stoull(s,&n);
    if (n!=s.size() || value==0 || value>std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("invalid orbital index: "+s);
    return static_cast<std::size_t>(value);
}
// Resolve an absolute MO index or a HOMO/LUMO-relative index (1-based).
std::size_t endpoint(const std::string& s,const Molecule& mol,Spin channel) {
    if (s.empty()) throw std::invalid_argument("empty orbital selector");
    if (std::isdigit(static_cast<unsigned char>(s[0]))) return number(s);
    std::size_t base=0, offset=0;
    if (s.rfind("homo",0)==0) {base=static_cast<std::size_t>(mol.occupied_count(channel));offset=4;}
    else if (s.rfind("lumo",0)==0) {base=static_cast<std::size_t>(mol.occupied_count(channel))+1;offset=4;}
    else throw std::invalid_argument("unknown MO selector: "+s);
    if (offset==s.size()) return base;
    const char sign=s[offset];
    if (sign!='+' && sign!='-') throw std::invalid_argument("invalid MO offset: "+s);
    const auto step=number(s.substr(offset+1));
    if (sign=='-' && step>=base) throw std::out_of_range("MO offset below first orbital: "+s);
    if (sign=='+') {
        if (step>std::numeric_limits<std::size_t>::max()-base) throw std::out_of_range("MO index overflow");
        return base+step;
    }
    return base-step;
}
// Expand inclusive ranges and de-duplicate MO selections in user order.
std::vector<std::size_t> select_mos(const std::string& selector,const Molecule& mol,Spin channel) {
    std::vector<std::size_t> result;
    std::set<std::size_t> seen;
    std::string expanded=selector=="frontier" ? "homo,lumo" : selector;
    std::size_t start=0;
    while (start<expanded.size()) {
        const auto comma=expanded.find(',',start);
        const auto token=expanded.substr(start,comma==std::string::npos ? std::string::npos : comma-start);
        const auto colon=token.find(':');
        auto first=endpoint(token.substr(0,colon),mol,channel);
        auto last=colon==std::string::npos ? first : endpoint(token.substr(colon+1),mol,channel);
        if (first>last || last>mol.orbital_count(channel))
            throw std::out_of_range("MO selector outside available orbitals: "+token);
        if (last-first>10000) throw std::invalid_argument("MO range too large");
        for (auto n=first;n<=last;++n) if (seen.insert(n).second) result.push_back(n);
        if (comma==std::string::npos) break;
        start=comma+1;
        if (start==expanded.size()) throw std::invalid_argument("trailing comma in MO selector");
    }
    if (result.empty()) throw std::invalid_argument("empty MO selector");
    return result;
}
std::string spin_name(Spin s) {return s==Spin::Alpha ? "alpha" : "beta";}
void list_mos(const Molecule& mol,Spin channel) {
    const auto count=mol.orbital_count(channel);
    const auto& energies=mol.orbital_energies(channel);
    const auto occupied=mol.occupied_count(channel);
    const bool restricted = !mol.unrestricted();
    std::cout << "Spin: " << (restricted ? "restricted" : spin_name(channel)) << "  orbitals: " << count << "  occupied: " << occupied << "\n";
    std::cout << " MO   Occupation   Energy (hartree)   Label\n";
    for (std::size_t i=1;i<=count;++i) {
        std::string label;
        if (i==static_cast<std::size_t>(occupied)) label="HOMO";
        if (i==static_cast<std::size_t>(occupied)+1) label="LUMO";
        std::cout << std::setw(4) << i << "   " << (i<=static_cast<std::size_t>(occupied) ? (restricted ? 2 : 1) : 0) << "            ";
        if (i<=energies.size()) std::cout << std::setw(15) << std::fixed << std::setprecision(8) << energies[i-1];
        else std::cout << std::setw(15) << "n/a";
        std::cout << "   " << label << "\n";
    }
}
void write_orbitals(const Molecule& mol,const chem::kumo::Grid& grid,const Options& o,Spin channel) {
    const auto indices=select_mos(o.selector,mol,channel);
    auto fields=mol.orbitals(grid,indices,channel);
    for (std::size_t j=0;j<indices.size();++j) {
        const auto n=indices[j];
        const auto filename=o.prefix+"_"+spin_name(channel)+"_mo_"+std::to_string(n)+".cube";
        fields[j].write_orbital_cube(filename,n,"Molecular orbital amplitude (signed), " + spin_name(channel) + " MO " + std::to_string(n));
        std::cout << filename << "\n";
    }
}
}
int main(int argc,char** argv) {
    try {
        if (argc==2 && std::string(argv[1])=="--help") {usage();return 0;}
        if (argc==2 && std::string(argv[1])=="--version") {std::cout << "kumo " CHEM_KUMO_VERSION "\n";return 0;}
        const auto o=parse(argc,argv);
        const auto mol=Molecule::from_fchk(o.input);
        if (o.list) {
            if (o.spin.empty() || o.spin=="alpha" || o.spin=="both") list_mos(mol,Spin::Alpha);
            if (o.spin=="beta" || o.spin=="both") list_mos(mol,Spin::Beta);
            if (!(o.spin.empty() || o.spin=="alpha" || o.spin=="beta" || o.spin=="both"))
                throw std::invalid_argument("--list-mo accepts alpha, beta, or both");
            return 0;
        }
        const auto grid=chem::kumo::Grid::around(mol.atoms(),o.margin,o.spacing);
        if (o.mo_mode) {
            if (!(o.spin.empty() || o.spin=="alpha" || o.spin=="beta" || o.spin=="both"))
                throw std::invalid_argument("--mo accepts --spin alpha, beta, or both");
            if (o.spin.empty() || o.spin=="alpha" || o.spin=="both") write_orbitals(mol,grid,o,Spin::Alpha);
            if (o.spin=="beta" || o.spin=="both") write_orbitals(mol,grid,o,Spin::Beta);
        } else {
            const auto spin=o.spin.empty() ? "total" : o.spin;
            std::vector<std::pair<Spin,std::string>> components;
            if (spin=="total" || spin=="all") components.emplace_back(Spin::Total,"density");
            if (spin=="alpha" || spin=="all") components.emplace_back(Spin::Alpha,"alpha");
            if (spin=="beta" || spin=="all") components.emplace_back(Spin::Beta,"beta");
            if (spin=="difference" || spin=="all") components.emplace_back(Spin::Difference,"spin");
            if (components.empty()) throw std::invalid_argument("density --spin must be total, alpha, beta, difference, or all");
            if (spin=="all") {
                auto fields=mol.spin_densities(grid);
                fields.total.write_cube(o.prefix+"_density.cube");
                fields.alpha.write_cube(o.prefix+"_alpha.cube");
                fields.beta.write_cube(o.prefix+"_beta.cube");
                fields.spin.write_cube(o.prefix+"_spin.cube");
            } else {
                mol.density(grid,components.front().first).write_cube(o.prefix+"_"+components.front().second+".cube");
            }
            for (const auto& [component,suffix]:components) {
                (void)component;
                std::cout << o.prefix << "_" << suffix << ".cube\n";
            }
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "kumo: " << e.what() << "\nTry 'kumo --help'.\n";
        return 1;
    }
}
