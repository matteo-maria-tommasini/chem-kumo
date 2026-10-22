# chem-kumo: MPFR s-primitive kernel audit

Input archive SHA-256: `4efe3f12ea90e0ed5a003b545b0debe89902872ad54414be5ca60beb8d10fff0`.

## Scope and result

Reviewed `local_mpfr_interval.hpp`, `local_mpfr_s_product.hpp`, `local_mpfr_s_pair_bound.hpp`, `local_mpfr_s_quartet_bound.hpp`, policy header and bundled tests. The tests reported by the developer on macOS are valuable evidence, but do not establish universal certification. **Do not mark the kernel `BoundReliability::Certified` yet.**

## Mathematical formula

For four isotropic s primitives with positive exponents, finite centers and finite coefficients, the absolute overlap equals

`|prod_i c_i| * (pi / sum_i alpha_i)^(3/2) * exp(-sum_{i<j} alpha_i alpha_j |R_i-R_j|^2 / sum_i alpha_i)`.

The pair-product formulation in the current quartet header is algebraically equivalent. Using lower interval endpoints for nonnegative decay terms, and an upper bound on the logarithmic volume and coefficient terms, yields an upper bound *provided every intermediate interval encloses the exact value*.

## Code findings

1. `sproduct_detail::mul` evaluates all four endpoint products with directed rounding; `div_positive` uses a reciprocal interval and rejects nonpositive lower denominators; `square` explicitly handles intervals crossing zero. These are mathematically appropriate for finite, ordered operands and nonexceptional MPFR operations.
2. `squartet_detail::product_data` builds Gaussian pair exponents and centers from exact binary64 input values through directed MPFR operations. It uses nonnegative distance-square terms and the correct lower bound of decay in the final negative exponent.
3. The quartet log-space expression computes the logarithmic volume with upper-rounded pi and lower-rounded total exponent. Coefficient logarithms are upper-rounded. The `denorm_min` threshold uses a downward-rounded logarithm; the `DBL_MAX` threshold uses an upward-rounded logarithm. These directions are conservative when preceding intermediates are valid.
4. **Critical unresolved exceptional-state gap:** MPFR operations within `product_data`, `mul`, `div_positive`, `square`, and the log-space construction are not systematically checked for NaN, overflow, underflow, or interval order. A directed lower bound of a positive quantity can underflow to zero; if subsequently used as a divisor the function throws. Other exceptional combinations may create NaN. A thrown exception is fail-closed only if every screening caller handles it by keeping the quartet.
5. **Separate API issue:** `s_primitive_product_intervals` replaces an upper amplitude of zero with `denorm_min` whenever both coefficients are nonzero. This is conservative if the exact amplitude is positive and no earlier invalid arithmetic occurred, but cannot repair an upstream invalid interval. Its behavior is not proof of general exception safety.
6. `s_primitive_quartet_overlap_upper` returns infinity for NaN log bounds, but may throw from `div_positive` on an exceptional interval. Its final `mpfr_get_d` result is not explicitly checked for NaN (although it checks the MPFR value). The API is not yet a total, fail-closed bound evaluator for all finite binary64 inputs.
7. MPFR's finite exponent range is vastly larger than binary64 but not unbounded. A 256-bit significand does not imply unlimited exponent range. Global exponent-range modifications are inappropriate in a library unless carefully scoped and restored; no such modifications are needed for the next milestone.
8. `local_bound_reliability.hpp` has an enum named `Certified`, but this does not itself supply a proof. Do not propagate that label from this kernel without an explicit reviewed contract.

## Closure criteria

- Make every numerical kernel return a typed outcome: `FiniteUpperBound`, `ConservativeInfinity`, or `InvalidInput`, and ensure that any exception/invalid interval leads to `Keep` in rigorous screening.
- Add invariant checks for `lo <= hi`, no NaN, and strictly positive lower bounds before division/logarithm; perform these checks at operation boundaries.
- Verify MPFR flag behavior in each stage with targeted injection tests and preserve callers' expected flag state where relevant.
- Confirm that exact-zero coefficients are handled before logarithms and that finite nonzero inputs never yield a zero upper bound.
- Run the existing MPFR tests, then full Meson suite on the Mac. Preserve `BoundReliability::Engineering` until code review of all paths is complete.

## Scope boundary

Even after closing the above items, certification applies only to absolute overlaps of **four s primitives** with exact binary64 input interpretation, not general electron-repulsion integrals, contracted shell quartets, Cartesian angular momenta, or density-weighted screening. Those require separate analytical bounds.
