# InstrumentOperator.jl
Collection of instrument line-shape methods for hyperspectral remote sensing. Basic objective is to provide composable routines to apply a variety of instrument line-shapes for the convolution and resampling of high resolution modeled radiance with an instrument operator. At the moment, we have FTS systems (sinc function convolved with assymetric box), spectrally varying lookup tables as for the Orbiting Carbon Observatory and generic distribution functions (e.g. Gaussian) as instrument kernels.

## Installation

InstrumentOperator can be installed using the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```julia
pkg> add https://github.com/RadiativeTransfer/InstrumentOperator.jl
```

## Example

Create an FTS instrument type with a maximum optical path difference (MOPD) of 2.5cm, a Field of View (FOV) of 7.9mrad and an assymetry factor of 0.02 (very close to GOSAT specifications:

```julia
# Create an instrument:
FTS = FTSInstrument(2.5, 1.9e-3, 0.02)
# Define a range for the ILS (in wavenumbers here)
x = -4:0.01:4
# Create a kernel at a center wavenumber of 13000cm⁻¹
FTSkernel = create_instrument_kernel(FTS, x,13000.0)
# plot
plot(x, FTSkernel.parent,label="FTS lineshape",lw=2)
```
![FTS_ils](https://user-images.githubusercontent.com/10467190/115968929-07cdd100-a4ef-11eb-839b-f22918743828.png)
