```@meta
EditURL = "<unknown>/docs/src/pages/tutorials/CrossSection_convolution.jl"
```

# Instrument Line Shapes

Using packages:

```@example CrossSection_convolution
using Plots
```

This needs to be installed from https://github.com/RadiativeTransfer/RadiativeTransfer.jl

```@example CrossSection_convolution
using RadiativeTransfer.Absorption
using InstrumentOperator
```

## Load HITRAN data and CO2 cross sections

```@example CrossSection_convolution
hitran_data = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=6000, ν_max=6400)
line_co2 = make_hitran_model(hitran_data, Voigt(), architecture=CPU())
```

## Compute a high resolution transmission spectrum

```@example CrossSection_convolution
Δν = 0.001
ν = 6275:Δν:6400;
nothing #hide
```

CO₂ cross section at 800hPa and 296K:

```@example CrossSection_convolution
σ_co2   = absorption_cross_section(line_co2, ν, 800.0, 296.0);
nothing #hide
```

Transmission through a typical atmosphere (8e21 molec/cm²)

```@example CrossSection_convolution
T = exp.(-8e21*σ_co2);
nothing #hide
```

plot high resolution transmission

```@example CrossSection_convolution
plot(ν, T, lw=2, label="High resolution")
```

Define the instrument kernel grid

```@example CrossSection_convolution
x = -8:Δν:8;
nothing #hide
```

## Create a kernel at a center wavenumber of 6300cm⁻¹
Use FTS with 5cm MOPD, FOV of 7.9mrad, no assymetry:

```@example CrossSection_convolution
FTS = FTSInstrument(5.0, 7.9e-3, 0.0)
FTSkernel = create_instrument_kernel(FTS, x,6300.0)
FTS_instr = FixedKernelInstrument(FTSkernel, collect(6280:0.01:6380))
```

## Convolve with instrument kernel

```@example CrossSection_convolution
T_conv = conv_spectra(FTS_instr, ν, T);
nothing #hide
```

overplot convolved transmission:

```@example CrossSection_convolution
plot!(FTS_instr.ν_out, T_conv, label="Convolved")
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

