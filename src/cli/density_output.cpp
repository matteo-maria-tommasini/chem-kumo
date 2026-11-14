#include "cli/density_output.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace kumo_cli {
void write_density_cubes(const Molecule& mol, const chem::kumo::Grid& grid,
                         const DensitySelection& selection, const Options& o, Console& console) {
    const auto& components = selection.components;
    const auto& properties = selection.properties;
    const double esp_tolerance = o.esp_tolerance;
    const auto& prefix = o.prefix;
    chem::kumo::Request request;
    for (const auto& [component, suffix] : components)
        for (const auto& property : properties)
            request.add(density_property(component, property));
    console.begin(evaluation_summary(components, properties));
    const auto results = chem::kumo::Evaluator(mol, evaluation_options(selection, o, console))
                             .evaluate(grid, request);
    console.end("Evaluated in");
    const auto files = components.size() * properties.size();
    console.begin("Writing " + std::to_string(files) + " Cube file" + (files > 1 ? "s" : ""));
    const auto write = [&](DensityComponent component, const std::string& suffix,
                           const DensityProperty& property) {
        auto values =
            scalar_values(results.get(density_property(component, property)), property.kind);
        const auto filename = prefix + "_" + suffix +
                              (property.kind == PropertyKind::Density ? "" : "_" + property.name) +
                              ".cube";
        const auto method = method_label(property, results.backend_used(), esp_tolerance);
        chem::kumo::Field(grid, mol.atoms(), std::move(values))
            .write_cube(filename,
                        property.description + ", " + component_label(component) + ", " + method);
        std::cout << filename << "\n";
    };
    for (const auto& property : properties)
        if (property.kind == PropertyKind::Density)
            for (const auto& [component, suffix] : components)
                write(component, suffix, property);
    for (const auto& [component, suffix] : components)
        for (const auto& property : properties)
            if (property.kind != PropertyKind::Density) write(component, suffix, property);
    console.end("Written in");
}
} // namespace kumo_cli
