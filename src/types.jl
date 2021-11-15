#####
##### Types for describing an Instrument Operations (convolution, resampling)
#####

"""
    type AbstractInstrumentOperator
Abstract AbstractInstrumentOperator type 
"""
abstract type AbstractInstrumentOperator end

"""
    type AbstractInstrument
Abstract Instrument type 
"""
abstract type AbstractInstrument end

"""
    struct FixedKernelInstrument{FT}

A struct which provides all parameters for the convolution with a kernel, which is identical across the instrument spectral axis

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FixedKernelInstrument{FT} <: AbstractInstrumentOperator
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
struct VariableKernelInstrument{FT,AA} <: AbstractInstrumentOperator
    "convolution Kernel" 
    kernel::OffsetArray{FT,2,AA}
    "Output spectral grid"
    ν_out::Array{FT,1}
    "Output index grid"
    ind_out::Array{Int,1}
end;

struct FTSInstrument{FT} <: AbstractInstrument
    "Maximum Optical Path Difference (cm)"
    MOPD::FT
    "Field of View (rad)"
    FOV::FT
    "Assymmetry parameter"
    β::FT
end

"""
    type AbstractNoiseModel
Abstract AbstractNoiseModel type 
"""
abstract type AbstractNoiseModel end

"""
    struct GratingNoiseModel{FT}

A struct which stores important variable to compute the noise of a grating spectrometer

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef  struct GratingNoiseModel <: AbstractNoiseModel
    "Integration time `[s]`" 
    t_int::Unitful.Time
    "Detector pixel size (assuming quadratic) `[μm]`"
    detector_size::Unitful.Length
    "Detector quantum efficiency"
    Qₑ
    "Efficiency of the optical bench (Telescope, grating and other optical elements)"
    η
    "F-number"
    Fnumber
    "Spectral Sampling Interval SSI or Δλ `[μm]`"
    Δλ::Unitful.Length
    "Readout noise `[e⁻]`"
    σ_read
    "Dark current `[e⁻/s]`"
    dark_current::PerTime
    # "Slit width `[μm]`"
    # slit_width::Unitful.Length
end;

function createGratingNoiseModel(ET, DS, FPA_qe, effTransmission, fnumber, SSI, RN, DC) 
    # @info "Creating instrument model"
    GratingNoiseModel(
       uconvert(u"s",ET),
       uconvert(u"μm",DS),
       FPA_qe, effTransmission,
       fnumber,
       uconvert(u"μm",SSI),
       RN,
       uconvert(u"s^-1",DC))
end

function Base.show(io::IO, m::GratingNoiseModel)
    compact = get(io, :compact, false)

    if !compact
        println("Instance of GratingNoiseModel:")
        println("Integration time         = ", m.t_int)
        println("Detector size            = ", m.detector_size)
        println("FPA Quantum efficiency   = ", m.Qₑ)
        println("Optical bench efficiency = ", m.η)
        println("F-number                 = ", m.Fnumber)
        println("Spectral Sampling        = ", m.Δλ)
        println("Readout noise            = ", m.σ_read)
        println("Dark_current             = ", m.dark_current)
    else
        show(io, "Instance of GratingNoiseModel")
    end
end

"""
    type AbstractL1File
Abstract AbstractL1File type 
"""
abstract type AbstractL1File end
struct L1_OCO <: AbstractL1File
    geometry::Dict{String, Any}
    ils::Dict{String, Any}
    meteo::Dict{String, Any}
    measurement::Dict{String, Any}
end;

abstract type AbstractMeasurement end

struct MeasurementOCO{FT} <: AbstractMeasurement
    "Wavelength axis"
    SpectralGrid::Array{FT,1}
    "Band ID axis"
    BandID::Array{Int,1}
    "Spectral radiances (concatenated)"
    SpectralMeasurement::Array{FT,1}
    "Radiance Unit"
    MeasurementUnit::String
    "Wavelength unit"
    GridUnit::String
    "Latitude [°]"
    latitude::FT
    "Longitude [°]"
    longitude::FT
    "Viewing Zenith Angle [°]"
    vza::FT
    "Solar Zenith Angle [°]"
    sza::FT
    "Relative Azimuth Angle [°]"
    raa::FT
    "Gravity at surface"
    g₀::FT
    "Polarization Angle [°]"
    ϕ::FT
    "Mueller Coefficients"
    mueller::Array{FT,1}
    "Pressure profile half levels"
    p::Array{FT,1}
    "Temperature profile full level"
    T::Array{FT,1}
    "Specific humidity full level"
    q::Array{FT,1}
    "Instrument Line Shape"
    ils::AbstractInstrumentOperator
end;

