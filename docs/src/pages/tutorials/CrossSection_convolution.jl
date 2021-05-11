# # Instrument Operator: Instrument Line Shapes

# ###  Using packages:
using Plots
using RadiativeTransfer.Absorption
using InstrumentOperator
using RadiativeTransfer.Architectures

# ## Load HITRAN data and CO2 cross sections
hitran_data = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=6000, ν_max=6400)
line_co2 = make_hitran_model(hitran_data, Voigt(), architecture=CPU())

# ## Compute a high resolution transmission spectrum 
Δν = 0.001
ν = 6275:Δν:6400;
# CO2 cross section at 800hPa and 296K:
σ_co2   = absorption_cross_section(line_co2, ν, 800.0, 296.0);
# Transmission through a typical atmosphere (8e21 molec/cm²)
T = exp.(-8e21*σ_co2);

# ### plot high resolution transmission
plot(ν, T, lw=2, label="High resolution")

# ## Define the instrument kernel grid
x = -8:Δν:8
# ## Create a kernel at a center wavenumber of 6300cm⁻¹
# Use FTS with 5cm MOPD, FOV of 7.9mrad, no assymetry:
FTS = FTSInstrument(5.0, 7.9e-3, 0.0)
FTSkernel = create_instrument_kernel(FTS, x,6300.0)
FTS_instr = FixedKernelInstrument(FTSkernel, collect(6280:0.01:6380))

# ## Convolve with instrument kernel
T_conv = conv_spectra(FTS_instr, ν, T);

# ### overplot convolved transmission:
plot!(FTS_instr.ν_out, T_conv, label="Convolved")