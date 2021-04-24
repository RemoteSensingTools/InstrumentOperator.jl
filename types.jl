#####
##### Types for describing an Instrument Operations (convolution, resampling)
#####

"""
    type InstrumentOperator
Abstract InstrumentOperator type 
"""
abstract type InstrumentOperator end

"""
    type Instrument
Abstract Instrument type 
"""
abstract type Instrument end

"""
    struct FixedKernelInstrument{FT}

A struct which provides all parameters for the convolution with a kernel, which is identical across the instrument spectral axis

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FixedKernelInstrument{FT} <: InstrumentOperator
    "convolution Kernel" 
    kernel::OffsetArray{FT,1}
    "Output spectral grid"
    ν_out::Array{FT,1}
end;

"""
    struct VariableInstrument{FT}

        A struct which provides all parameters for the convolution with a kernel, which varies across the instrument spectral axis

# Fields
$(DocStringExtensions.FIELDS)
"""
struct VariableKernelInstrument{FT,AA} <: InstrumentOperator
    "convolution Kernel" 
    kernel::OffsetArray{FT,2,AA}
    "Output spectral grid"
    ν_out::Array{FT,1}
    "Output index grid"
    ind_out::Array{Int,1}
end;

struct FTSInstrument{FT} <: Instrument
    "Maximum Optical Path Difference (cm)"
    MOPD::FT
    "Field of View (rad)"
    FOV::FT
    "Assymmetry parameter"
    β::FT
end
