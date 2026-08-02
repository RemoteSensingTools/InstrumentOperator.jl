# Changelog

## v0.1.3 - 2026-08-01

### Added

- `CompactVariableKernelInstrument`, a grid-bound, compact representation of
  spectrally varying tabulated instrument line shapes.
- `prepare_compact_ils` for conservative trapezoid quadrature, per-channel
  normalization, and compact contiguous support construction.
- Allocation-free `conv_spectra!` methods for individual spectra and batched
  design matrices, plus allocating `conv_spectra` wrappers.
- Aqua and the compatible JET generation as release checks on Julia 1.10,
  1.11, and 1.12.

### Changed

- Added `LinearAlgebra` as an explicit dependency for optimized compact dot
  products.
- Removed the stale `UnitfulRecipes` dependency.
- Invalid or missing ILS table inputs now raise an exception instead of
  returning ambiguous `nothing` tuples.

### Validation

- Compact weights preserve a spectrally flat input exactly to floating-point
  precision.
- OCO-2 tests agree with independently recomputed exact quadrature to below
  one part per million while reducing a hot convolution to 5--7 microseconds
  with zero allocations.
