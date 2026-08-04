# Changelog

## v0.1.4 - 2026-08-03

### Added

- Optional CUDA support for `CompactVariableKernelInstrument`, loaded only
  when CUDA.jl is present.
- `gpu_operator` for moving a prepared compact operator to an NVIDIA GPU.
- Separate CUDA correctness tests and a batched OCO-style benchmark.

### Changed

- Compact vector and matrix convolution kernels are implemented with
  KernelAbstractions.jl. The numerical kernels are backend-independent; the
  CUDA extension only provides device storage and backend selection.
- Compact operator array fields now retain their concrete host or device array
  types, and allocating convolution preserves the input array backend.

### Validation

- CUDA and CPU results agree within `2e-6` in the dedicated Float32 test.
- On an NVIDIA L40S, an OCO-style 1,551-point to 101-channel convolution
  sustains approximately 64 million spectra per second for an 8,192-spectrum
  batch.

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
