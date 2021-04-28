[![dev][docs-latest-img]][docs-latest-url] 

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://radiativetransfer.github.io/InstrumentOperator.jl/


# InstrumentOperator.jl
Collection of instrument line-shape methods for hyperspectral remote sensing. Basic objective is to provide composable routines to apply a variety of instrument line-shapes (ILS) for the convolution and resampling of high resolution modeled radiance with an instrument operator. At the moment, we have FTS systems (sinc function convolved with assymetric box), spectrally varying lookup tables as for the Orbiting Carbon Observatory and generic distribution functions (e.g. Gaussian) as instrument kernels. 


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

We can also use Julia's `Distributions` package and create custom kernels with Continous distributions and convolutions thereof, e.g. an AVIRIS-NG like ILS (5nm pixels):
```julia
# Define a range for the ILS (resolution always has to be the same as your high resolution spectrum you want to convolve!).
x = -15:0.01:15
# Create a custom kernel using two Distributions, convolution of box and Gaussian (e.g. pixel width as box, smoothed by Gaussian)
avNGkernel = create_instrument_kernel([Normal(0, 1),Uniform(-2.5,2.5)], x)
```
![AV_ils](https://user-images.githubusercontent.com/10467190/116157489-ec56f780-a6a1-11eb-9048-8b6938fec69f.png)



