// Gzip-compressed formatted checkpoint (FCHK) files: for O2 and acrolein, the
// compressed file and a decompressed copy give identical scalars, coordinates
// and density derivatives (exact equality); compression is detected from the
// file content, not the suffix; a damaged or truncated gzip file is rejected.
// Without a "Nuclear charges" record, the charges fall back to the atomic
// numbers, and an effective core potential file without it is rejected. For
// the three hexabenzocoronene (HBC) files, the atom, electron and basis
// function counts are checked.
//
// Usage: kumo-fchk-gzip-test O2.fchk.gz ACROLEIN.fchk.gz
//            HBC1.fchk.gz HBC2.fchk.gz HBC3.fchk.gz
#include <chem/kumo.hpp>
#include <chem/kumo/fchk.hpp>
#include "molecule_test_access.hpp"
#include "test_support.hpp"

#include <zlib.h>

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

// Materialize a plain-text copy only in the temporary test directory.
// The committed fixture remains gzip-compressed.
void write_plain(const std::filesystem::path& input, const std::filesystem::path& output) {
    gzFile source = gzopen(input.string().c_str(), "rb");
    require(source != nullptr, "cannot open reference gzip");
    std::ofstream target(output, std::ios::binary);
    require(static_cast<bool>(target), "cannot create temporary plain FCHK");
    std::array<char, 32768> bytes{};
    int count = 0;
    while ((count = gzread(source, bytes.data(), static_cast<unsigned>(bytes.size()))) > 0) {
        target.write(bytes.data(), count);
        require(static_cast<bool>(target), "cannot write temporary plain FCHK");
    }
    require(count == 0 && gzclose(source) == Z_OK, "cannot decompress reference gzip");
    target.close();
    require(static_cast<bool>(target), "cannot close temporary plain FCHK");
}

void expect_failure(const std::filesystem::path& path) {
    try {
        (void)chem::kumo::FormattedCheckpoint(path);
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error("corrupt gzip was accepted: " + path.string());
}

// Without a "Nuclear charges" record the nuclear charges fall back to the
// atomic numbers, which is right for all-electron files. With an ECP record the
// effective charges are unknown, so loading must fail.
void verify_nuclear_charge_fallback(const std::filesystem::path& plain,
                                    const std::filesystem::path& directory) {
    std::ifstream input(plain, std::ios::binary);
    const std::string text((std::istreambuf_iterator<char>(input)),
                           std::istreambuf_iterator<char>());
    // Drop the "Nuclear charges" header and its data lines (which start with a
    // space, while every record header starts with a letter).
    const auto begin = text.find("\nNuclear charges ");
    require(begin != std::string::npos, "fixture has no Nuclear charges record");
    auto end = text.find('\n', begin + 1);
    while (end != std::string::npos && end + 1 < text.size() && text[end + 1] == ' ')
        end = text.find('\n', end + 1);
    require(end != std::string::npos, "Nuclear charges record is not terminated");
    const std::string without_charges = text.substr(0, begin) + text.substr(end);

    const auto all_electron = directory / (plain.filename().string() + ".no-charges.fchk");
    std::ofstream(all_electron, std::ios::binary) << without_charges;
    const auto mol = chem::kumo::Molecule::from_fchk(all_electron);
    for (const auto& atom : mol.atoms())
        kumo_test::check(atom.nuclear_charge == static_cast<double>(atom.atomic_number),
                         plain.filename().string() +
                             ": without Nuclear charges, the nuclear charge is Z = " +
                             std::to_string(atom.atomic_number) + " (got " +
                             kumo_test::num(atom.nuclear_charge) + ")");

    const auto with_ecp = directory / (plain.filename().string() + ".ecp.fchk");
    std::ofstream(with_ecp, std::ios::binary)
        << without_charges << "ECP-LMax                                   I   N=           1\n"
        << "           0\n";
    kumo_test::check_throws<std::runtime_error>(
        [&] { (void)chem::kumo::Molecule::from_fchk(with_ecp); },
        plain.filename().string() + ": an ECP file without Nuclear charges is rejected");
}

void verify(const std::filesystem::path& zipped, const std::filesystem::path& directory) {
    require(zipped.extension() == ".gz", "expected a compressed FCHK fixture");
    const auto plain_path = directory / zipped.stem();
    write_plain(zipped, plain_path);
    const auto plain = chem::kumo::FormattedCheckpoint(plain_path);
    const auto compressed = chem::kumo::FormattedCheckpoint(zipped);
    for (const auto* label : {"Number of atoms", "Number of basis functions",
                              "Number of alpha electrons", "Number of beta electrons"}) {
        require(plain.integer(label) == compressed.integer(label), "FCHK scalar mismatch");
    }
    const auto& a = plain.reals("Current cartesian coordinates");
    const auto& b = compressed.reals("Current cartesian coordinates");
    require(a == b, "FCHK coordinates mismatch");
    const auto mol_a = chem::kumo::Molecule::from_fchk(plain_path);
    const auto mol_b = chem::kumo::Molecule::from_fchk(zipped);
    verify_nuclear_charge_fallback(plain_path, directory);
    const Eigen::Vector3d point(0.2, -0.3, 0.4);
    for (const auto spin :
         {chem::kumo::DensityComponent::Total, chem::kumo::DensityComponent::Alpha,
          chem::kumo::DensityComponent::Beta, chem::kumo::DensityComponent::Difference}) {
        const auto x = kumo_test::Kernels::density_derivatives(mol_a, point, spin);
        const auto y = kumo_test::Kernels::density_derivatives(mol_b, point, spin);
        require(x.value == y.value && x.laplacian == y.laplacian && x.tau == y.tau &&
                    (x.gradient - y.gradient).norm() == 0.0 &&
                    (x.hessian - y.hessian).norm() == 0.0,
                "compressed density derivatives mismatch");
    }
    // Compression detection uses magic bytes rather than filename suffix.
    const auto renamed = directory / (zipped.stem().string() + ".dat");
    std::filesystem::copy_file(zipped, renamed, std::filesystem::copy_options::overwrite_existing);
    require(chem::kumo::FormattedCheckpoint(renamed).integer("Number of atoms") ==
                plain.integer("Number of atoms"),
            "gzip magic detection failed");

    // Deliberately damage the gzip trailer (CRC32/ISIZE) and truncate the file.
    const auto corrupt = directory / (zipped.stem().string() + ".bad.gz");
    std::filesystem::copy_file(zipped, corrupt, std::filesystem::copy_options::overwrite_existing);
    const auto size = std::filesystem::file_size(corrupt);
    require(size > 8, "unexpected gzip fixture size");
    {
        std::fstream stream(corrupt, std::ios::binary | std::ios::in | std::ios::out);
        stream.seekg(static_cast<std::streamoff>(size - 8));
        char original = 0;
        stream.read(&original, 1);
        require(static_cast<bool>(stream), "cannot read gzip trailer");
        stream.seekp(static_cast<std::streamoff>(size - 8));
        const char changed = static_cast<char>(static_cast<unsigned char>(original) ^ 0x80u);
        stream.write(&changed, 1);
        require(static_cast<bool>(stream), "cannot modify gzip trailer");
    }
    expect_failure(corrupt);
    std::filesystem::resize_file(corrupt, size - 5);
    expect_failure(corrupt);
    std::cout << zipped.filename().string() << " gzip/plain equality and corruption checks OK\n";
}

void verify_large_reference(const std::filesystem::path& zipped) {
    const auto fchk = chem::kumo::FormattedCheckpoint(zipped);
    require(fchk.integer("Number of atoms") == 60, "HBC atom count mismatch");
    require(fchk.integer("Number of alpha electrons") == 135, "HBC alpha count mismatch");
    require(fchk.integer("Number of beta electrons") == 135, "HBC beta count mismatch");
    const auto name = zipped.filename().string();
    const int expected_basis = name.find("def2tzvp") != std::string::npos  ? 1410
                               : name.find("def2svp") != std::string::npos ? 678
                                                                           : 720;
    require(fchk.integer("Number of basis functions") == expected_basis,
            "HBC basis count mismatch");
    std::cout << name << " compressed reference OK\n";
}

} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(
            argc, argv, 5, 5, "O2.fchk.gz ACROLEIN.fchk.gz HBC1.fchk.gz HBC2.fchk.gz HBC3.fchk.gz"))
        return 2;
    const auto directory =
        std::filesystem::temp_directory_path() /
        ("kumo-gzip-test-" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory);
    try {
        verify(argv[1], directory);
        verify(argv[2], directory);
        for (int i = 3; i < argc; ++i)
            verify_large_reference(argv[i]);
        std::filesystem::remove_all(directory);
        return kumo_test::failures == 0 ? 0 : 1;
    } catch (const std::exception& e) {
        std::filesystem::remove_all(directory);
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
