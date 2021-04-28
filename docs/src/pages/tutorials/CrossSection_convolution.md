```@meta
EditURL = "<unknown>/docs/src/pages/tutorials/CrossSection_convolution.jl"
```

```@example CrossSection_convolution
using Pkg
using Plots
```

Add packages:

```@example CrossSection_convolution
Pkg.add(url="https://github.com/RadiativeTransfer/Architectures.jl")
Pkg.add(url="https://github.com/RadiativeTransfer/Absorption.jl")
push!(LOAD_PATH,"./src/")

using Absorption
using InstrumentOperator
using Architectures

hitran_data = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=6000, ν_max=6400)
line_co2 = make_hitran_model(hitran_data, Voigt(), architecture=CPU())

Δν = 0.001
ν = 6275:Δν:6400;
σ_co2   = absorption_cross_section(line_co2, ν, 200.0     , 296.0);
nothing #hide
```

Transmission through a typical atmosphere

```@example CrossSection_convolution
T = exp.(-8e21*σ_co2)

plot(ν, T, lw=2)


x = -8:Δν:8
```

Create a kernel at a center wavenumber of 6300cm⁻¹

```@example CrossSection_convolution
FTS = FTSInstrument(45.0, 1.2e-3, 0.0)
FTSkernel = create_instrument_kernel(FTS, x,6300.0)
FTS_instr = FixedKernelInstrument(FTSkernel, collect(6280:0.01:6380))

T_conv = conv_spectra(FTS_instr, ν, T);
plot!(FTS_instr.ν_out, T_conv)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

