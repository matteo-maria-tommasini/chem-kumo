#pragma once

namespace chem::kumo {

/// Component of the electron density: the total density, the alpha or beta
/// spin density, or the spin-difference density alpha - beta (which can be
/// negative).
enum class DensityComponent { Total, Alpha, Beta, Difference };

/// Spin channel of a molecular orbital (MO). A restricted calculation has one
/// set of spatial MOs, used for both channels.
enum class SpinChannel { Alpha, Beta };

} // namespace chem::kumo
