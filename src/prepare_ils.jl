
"""
$(FUNCTIONNAME)(grid_x::AbstractRange, ils_in::Array{FT}, ils_Δ::Array{FT}, extended_dims::Array{Int}=[]) where {FT <: AbstractFloat}

Pre-compute the ILS table input as function of spectral distance from center converted to the modeling grid 
Input: grid_x, ils_in, ils_Δ, extended_dims
Output: Offset-Array with tabulated responses interpolated to grid_x
"""
function prepare_ils_table(
    grid_x::AbstractRange,
    ils_in::Array{FT},
    ils_Δ::Array{FT},
    extended_dims::Array{Int}=Int[],
) where {FT <: AbstractFloat}
    @assert minimum(abs.(grid_x)) < eps(FT) "grid_x must include 0 (center pixel)"
    ind_0 = argmin(abs.(grid_x))
    axis_pixel = (-ind_0 + 1):(grid_x.len - ind_0)
    # Dimension of ILS table (at least 1D, first dimension needs to be across wavenumber/wavelength)
    dims  = size(ils_in);
    # number of spectral positions of ILS table
    n_x = dims[1]
    # Number of ILS per detector position (can be 1 if ILS is constant across detector grid)
    n_pos = dims[2]
    ils    = view(ils_in, :, :,  extended_dims...);
    ils_Δ_ = view(ils_Δ,  :, :,  extended_dims...);
    ils_pixel = zeros(FT, grid_x.len, n_pos);

    for i = 1:n_pos
        ind = findall(minimum(ils_Δ_[:,i]) .< grid_x .< maximum(ils_Δ_[:,i]));
        interp = Interpolations.linear_interpolation(ils_Δ_[:, i], ils[:, i])
        ils_pixel[ind,i] = interp.(grid_x[ind]);
    end
    # normalize here
    #return OffsetArray(ils_pixel ./ sum(ils_pixel, dims=1), axis_pixel, 1:n_pos)
    return OffsetArray(ils_pixel * FT(grid_x.step) , axis_pixel, 1:n_pos)
end

@inline function linear_table_value(
    coordinate::AbstractVector,
    value::AbstractVector,
    query,
)
    upper = searchsortedfirst(coordinate, query)
    upper = clamp(upper, firstindex(coordinate) + 1, lastindex(coordinate))
    lower = upper - 1
    fraction = (query - coordinate[lower]) /
        (coordinate[upper] - coordinate[lower])
    return value[lower] + fraction * (value[upper] - value[lower])
end

function compact_half_widths(
    ils_delta::AbstractMatrix{FT},
    half_width,
) where {FT<:AbstractFloat}
    number_of_channels = size(ils_delta, 2)
    if half_width === nothing
        widths = Vector{FT}(undef, number_of_channels)
        for channel in 1:number_of_channels
            delta = view(ils_delta, :, channel)
            widths[channel] = max(abs(first(delta)), abs(last(delta)))
        end
        return widths
    elseif half_width isa Real
        return fill(FT(half_width), number_of_channels)
    elseif half_width isa AbstractVector
        length(half_width) == number_of_channels || throw(DimensionMismatch(
            "half_width must have one value per output channel",
        ))
        return FT.(half_width)
    end
    throw(ArgumentError("half_width must be nothing, a scalar, or a vector"))
end

"""
    prepare_compact_ils(
        ν_in,
        ν_out,
        ils_response,
        ils_delta;
        half_width=nothing,
    )

Prepare an exact, compact, conservative variable-ILS operator on the fixed
high-resolution grid `ν_in`.

`ils_response` and `ils_delta` have dimensions `(ILS samples, output
channels)`. For each detector center in `ν_out`, the tabulated response is
linearly interpolated onto `ν_in`, multiplied by trapezoid cell widths, and
normalized. Only its contiguous support is retained. The resulting
[`CompactVariableKernelInstrument`](@ref) can apply one spectrum or a matrix
whose columns are spectra.

By default, each channel uses the larger absolute endpoint of its own ILS
offset table as its half-width. A scalar `half_width` reproduces an instrument
with one common integration half-width; a vector supplies one value per output
channel. Linear extrapolation at the first grid point just outside a requested
bound matches the ReFRACtor grating-convolution convention.

The prepared operator is valid only for this exact input grid. Both grids and
every ILS offset table must be finite and strictly increasing.
"""
function prepare_compact_ils(
    ν_in::AbstractVector{FT},
    ν_out::AbstractVector{FT},
    ils_response::AbstractMatrix{FT},
    ils_delta::AbstractMatrix{FT};
    half_width=nothing,
) where {FT<:AbstractFloat}
    axes(ν_in, 1) == Base.OneTo(length(ν_in)) || throw(ArgumentError(
        "the high-resolution input grid must use one-based indexing",
    ))
    axes(ν_out, 1) == Base.OneTo(length(ν_out)) || throw(ArgumentError(
        "the output grid must use one-based indexing",
    ))
    size(ils_response) == size(ils_delta) || throw(DimensionMismatch(
        "ILS response and offset arrays must have identical dimensions",
    ))
    size(ils_response, 1) >= 2 || throw(ArgumentError(
        "each ILS table requires at least two samples",
    ))
    size(ils_response, 2) == length(ν_out) || throw(DimensionMismatch(
        "ILS tables must have one column per output channel",
    ))
    length(ν_in) >= 2 || throw(ArgumentError(
        "the high-resolution input grid requires at least two samples",
    ))
    isempty(ν_out) && throw(ArgumentError("the output grid cannot be empty"))
    all(isfinite, ν_in) || throw(ArgumentError(
        "the high-resolution input grid must be finite",
    ))
    all(isfinite, ν_out) || throw(ArgumentError(
        "the output grid must be finite",
    ))
    all(isfinite, ils_delta) || throw(ArgumentError(
        "ILS offsets must be finite",
    ))
    all(isfinite, ils_response) || throw(ArgumentError(
        "ILS responses must be finite",
    ))
    all(diff(ν_in) .> zero(FT)) || throw(ArgumentError(
        "the high-resolution input grid must be strictly increasing",
    ))
    all(diff(ν_out) .> zero(FT)) || throw(ArgumentError(
        "the output grid must be strictly increasing",
    ))

    number_of_channels = length(ν_out)
    for channel in 1:number_of_channels
        all(diff(view(ils_delta, :, channel)) .> zero(FT)) ||
            throw(ArgumentError(
                "ILS offsets must be strictly increasing in channel $(channel)",
            ))
    end
    widths = compact_half_widths(ils_delta, half_width)
    all(isfinite, widths) && all(widths .> zero(FT)) || throw(ArgumentError(
        "ILS half-widths must be finite and positive",
    ))

    first_indices = Vector{Int}(undef, number_of_channels)
    lengths = Vector{Int}(undef, number_of_channels)
    maximum_length = 0
    for channel in 1:number_of_channels
        center = ν_out[channel]
        first_index = searchsortedfirst(ν_in, center - widths[channel])
        last_index = searchsortedfirst(ν_in, center + widths[channel])
        1 <= first_index < last_index <= length(ν_in) || throw(ArgumentError(
            "input grid does not contain the full ILS for output channel $(channel)",
        ))
        first_indices[channel] = first_index
        lengths[channel] = last_index - first_index + 1
        maximum_length = max(maximum_length, lengths[channel])
    end

    weights = zeros(FT, maximum_length, number_of_channels)
    for channel in 1:number_of_channels
        count = lengths[channel]
        first_index = first_indices[channel]
        center = ν_out[channel]
        delta = view(ils_delta, :, channel)
        response = view(ils_response, :, channel)
        local_grid = view(ν_in, first_index:(first_index + count - 1))
        for local_index in 1:count
            weights[local_index, channel] = linear_table_value(
                delta,
                response,
                local_grid[local_index] - center,
            )
        end
        weights[1, channel] *= (local_grid[2] - local_grid[1]) / FT(2)
        for local_index in 2:(count - 1)
            weights[local_index, channel] *= (
                local_grid[local_index + 1] -
                local_grid[local_index - 1]
            ) / FT(2)
        end
        weights[count, channel] *=
            (local_grid[count] - local_grid[count - 1]) / FT(2)
        normalization = sum(view(weights, 1:count, channel))
        isfinite(normalization) && normalization > zero(FT) ||
            throw(ArgumentError(
                "ILS response has non-positive integral in channel $(channel)",
            ))
        weights[1:count, channel] ./= normalization
    end

    return CompactVariableKernelInstrument(
        weights,
        first_indices,
        lengths,
        ν_in isa Vector{FT} ? ν_in : collect(ν_in),
        ν_out isa Vector{FT} ? ν_out : collect(ν_out),
    )
end
# This can still derive min/max range automatically, need to double check

function create_instrument_kernel(FTS::FTSInstrument, grid_x::AbstractRange, ν̅)
    @unpack MOPD, β, FOV = FTS
    # Need to make sure this is 0!
    i₀ = argmin(abs.(grid_x))
    axis_pixel = (-i₀ + 1):(grid_x.len - i₀)
    @info "Center of OffsetArrays is " minimum(abs.(grid_x))
    #Δσ = 1/(2FTS.OPD)
    
    # Create sinc FTS kernel
    XX = 2MOPD * grid_x
    sinc_kernel = OffsetArray(2MOPD * sinc.(XX),axis_pixel)
    
    # Now to FOV box:
    Θ = ν̅/2 * FOV^2
    β = β * Θ
    @show Θ, β
    # Create empty box_kernel:
    box = collect(0*grid_x);
    for i in eachindex(grid_x)
        if grid_x[i] ≤ 0
            if Θ > β && grid_x[i] ≥ -(Θ - β)
                box[i] = 1
            elseif grid_x[i] ≥ -(Θ + β)
                box[i] = acos((grid_x[i]^2 + β^2 - Θ^2)/(2grid_x[i] * β)) / π
            end
        end
    end
    box_kernel = OffsetArray(box,axis_pixel)
    box_kernel /= sum(box_kernel)
    sinc_kernel /= sum(sinc_kernel)
    return imfilter(box_kernel,sinc_kernel)
end

"Create kernel from a single continous distribution function"
function create_instrument_kernel(di::ContinuousUnivariateDistribution, grid_x::AbstractRange)
    ils = pdf.(di,grid_x)
    i₀ = argmin(abs.(grid_x))
    axis_pixel = (-i₀ + 1):(grid_x.len - i₀)
    ils /= sum(ils)
    return OffsetArray(ils,axis_pixel)
end

"Create kernel from multiple single continous distribution function (convolves those)"
function create_instrument_kernel(di::Array{T}, grid_x::AbstractRange) where T <: ContinuousUnivariateDistribution
    ils = pdf.(di[1],grid_x)
    for i=2:length(di)
        # might change this in future to just do multiplications in fft space and then ifft back
        ils = imfilter(ils, pdf.(di[i],grid_x))
    end
    i₀ = argmin(abs.(grid_x))
    axis_pixel = (-i₀ + 1):(grid_x.len - i₀)
    ils /= sum(ils)
    return OffsetArray(ils,axis_pixel)
end
