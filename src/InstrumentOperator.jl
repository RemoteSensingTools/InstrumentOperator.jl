module InstrumentOperator

using OffsetArrays
using Polynomials
using NCDatasets
using Interpolations
using Parameters
using DocStringExtensions
using Distributions
using ImageFiltering   # for convolutions
using Unitful
using OrderedCollections
using UnitfulEquivalences
using YAML

@derived_dimension PerTime inv(Unitful.𝐓)

include("types.jl")                   # All types used in this module
include("ils_tables_io.jl")           # IO functions
include("prepare_ils.jl")             # function for ILS preparation
include("instrument_convolutions.jl") # conv using ImageFiltering, can be changed later
include("noise_functions.jl")
include("io_L1_files.jl")

export FixedKernelInstrument, VariableKernelInstrument, FTSInstrument, create_instrument_kernel, conv_spectra   

end
