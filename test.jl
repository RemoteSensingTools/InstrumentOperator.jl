using OffsetArrays
using Polynomials
using NCDatasets
using Interpolations
using DocStringExtensions
using Plots
using Parameters
# using DataInterpolations



include("types.jl")
include("ils_tables_io.jl")
include("prepare_ils.jl")
include("instrument_convolutions.jl")

FT = Float32

using NCDatasets
using Polynomials
file = "/net/fluo/data2/groupMembers/cfranken/oco2_L1bScGL_18689a_180105_B10003r_200316045235.h5"
ils_Δ, ils_in, dispersion = read_ils_table(file, "src/InstrumentOperators/json/ils_oco2.json");
band = 1
footprint = 1
# Footprint and band index:
extended_dims  = [footprint,band]
dispPoly = Polynomial(view(dispersion, :, extended_dims...))

ν = FT.(dispPoly.(0:1015))
# resolution (μm)
res = FT(0.001/1e3)
# Max range of ILS from center (µm)
kernel_range = FT(0.5/1e3)

# Kernel grid
grid_x = -kernel_range:res:kernel_range
ind_out = collect(1:1016);

ils_pixel   = prepare_ils_table(grid_x, ils_in, ils_Δ,extended_dims)
oco2_kernel = VariableKernelInstrument(ils_pixel, ν, ind_out)

ν_mod = collect(FT, 0.758:Float32(res / 1e3):0.771)
spectrum = zeros(FT, length(ν_mod))
spectrum[length(spectrum) ÷ 2]=1;
y_conv =  conv_spectra(oco2_kernel, ν_mod, spectrum);