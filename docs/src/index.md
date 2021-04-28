# InstrumentOperator.jl
Collection of instrument line-shape methods for hyperspectral remote sensing. Basic objective is to provide composable routines to apply a variety of instrument line-shapes (ILS) for the convolution and resampling of high resolution modeled radiance (R) with an instrument operator. At the moment, we have FTS systems (sinc function convolved with assymetric box), spectrally varying lookup tables as for the Orbiting Carbon Observatory and generic distribution functions (e.g. Gaussian) as instrument kernels. 

The main purpose is thus to provide a flexible framework to create ILS functions as ``ILS(\Delta\lambda)`` with appropriate upper/lower limits to perform the convolution at the high resolution grid:
```math
(R*ILS)(\lambda_{hr}) = \int_{-\infty}^{\infty} R(x)ILS(\lambda_{hr}-x)dx
```
and the necessary interpolation from the high resolution (hr) grid to the low resolution (lr) detector grid after convolution. 
```math
(R*ILS)(\lambda_{hr}) \xrightarrow{\text{Cubic Spline}} (R*ILS)(\lambda_{lr})
```
To simplify calculations, we require the high resolution grid to be equidistant, ideally defined as Range.  

## Installation

InstrumentOperator can be installed using the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```julia
pkg> add https://github.com/RadiativeTransfer/InstrumentOperator.jl
```
