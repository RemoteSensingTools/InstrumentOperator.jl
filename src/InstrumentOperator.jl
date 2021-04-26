module InstrumentOperator
using OffsetArrays
using Polynomials
using NCDatasets
using Interpolations
using Parameters
using DocStringExtensions
# using DataInterpolations

include("types.jl")                # All types used in this module
include("ils_tables_io.jl")        # IO functions
include("prepare_ils.jl")          # function for ILS preparation

end
