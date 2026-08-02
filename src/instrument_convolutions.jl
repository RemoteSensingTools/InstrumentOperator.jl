
"Convolves and resamples the input spectrum with a fixed kernel"
function conv_spectra(m::FixedKernelInstrument, ν, spectrum)
    s = imfilter(spectrum, m.kernel)
    interp_cubic = cubic_spline_interpolation(spectral_range(ν), s; extrapolation_bc=Interpolations.Flat())
    return interp_cubic(m.ν_out)
end;

"Convolves and resamples the input spectrum with a variable kernel (per spectral position)"
function conv_spectra(m::VariableKernelInstrument, ν, spectrum; stride=1)
    FT = eltype(m.ν_out)
    FT2 = eltype(spectrum)
    # Define grid where to perform convolution:
    
    # Padding at both sides required:
    start = argmin(abs.(ν.-m.ν_out[1]))-1
    stop  = argmin(abs.(ν.-m.ν_out[end]))+1
    # @show start, stop
    off = ceil(Int, size(m.kernel, 1) / 2)
    @assert off ≤ start "Start range of model grid too close to output grid (needs buffer) $(off), $(start)"
    @assert off ≤ length(ν)-stop "Stop range of model grid too close to output grid (needs buffer) $(off), $(length(ν)-stop)"
    ind = start:stride:stop
    
    # knots where convolution will be applied to
    knots = view(ν, ind)
    te = linear_interpolation(spectral_range(m.ν_out), FT.(m.ind_out); extrapolation_bc=Interpolations.Flat())
    spec_out = zeros(FT2, length(knots));
    for i in eachindex(knots)
        # Simple first, nearest neighbor ILS
        ind_fraction = round(Int, te(knots[i]));
        kernel = view(m.kernel, :, ind_fraction)
        for j in eachindex(kernel)
            spec_out[i] += kernel[j] * spectrum[ind[i] + j] 
        end
    end
    # Change this later to only perform conv around output grid!
    fin = linear_interpolation(spectral_range(ν[ind]), spec_out; extrapolation_bc=Interpolations.Flat())
    return fin(m.ν_out)
end;

function check_compact_grid(m::CompactVariableKernelInstrument, ν)
    axes(ν, 1) == Base.OneTo(length(ν)) || throw(ArgumentError(
        "the high-resolution input grid must use one-based indexing",
    ))
    length(ν) == length(m.ν_in) || throw(DimensionMismatch(
        "input grid length does not match the prepared compact ILS operator",
    ))
    if ν !== m.ν_in
        @inbounds for index in eachindex(ν, m.ν_in)
            isequal(ν[index], m.ν_in[index]) || throw(ArgumentError(
                "input grid does not match the grid used to prepare the compact ILS operator",
            ))
        end
    end
    return nothing
end

"""
    conv_spectra!(output, instrument::CompactVariableKernelInstrument, ν, spectrum)

Apply a prepared compact variable ILS without allocating. `spectrum` can be a
vector, or a matrix whose columns are independent spectra. The corresponding
`output` shape is `(length(instrument.ν_out),)` or
`(length(instrument.ν_out), size(spectrum, 2))`.

The input grid must exactly match the grid used by [`prepare_compact_ils`](@ref).
Input and output arrays must not alias.
"""
function conv_spectra!(
    output::AbstractVector,
    m::CompactVariableKernelInstrument,
    ν::AbstractVector,
    spectrum::AbstractVector,
)
    check_compact_grid(m, ν)
    axes(spectrum, 1) == Base.OneTo(length(spectrum)) || throw(ArgumentError(
        "the input spectrum must use one-based indexing",
    ))
    axes(output, 1) == Base.OneTo(length(output)) || throw(ArgumentError(
        "the output spectrum must use one-based indexing",
    ))
    length(spectrum) == length(ν) || throw(DimensionMismatch(
        "spectrum and input grid lengths do not match",
    ))
    length(output) == length(m.ν_out) || throw(DimensionMismatch(
        "output length does not match the compact ILS detector grid",
    ))
    Base.mightalias(output, spectrum) && throw(ArgumentError(
        "compact ILS input and output arrays must not alias",
    ))
    @inbounds for channel in eachindex(output)
        first_index = m.first_indices[channel]
        count = m.lengths[channel]
        output[channel] = dot(
            view(m.weights, 1:count, channel),
            view(spectrum, first_index:(first_index + count - 1)),
        )
    end
    return output
end

function conv_spectra!(
    output::AbstractMatrix,
    m::CompactVariableKernelInstrument,
    ν::AbstractVector,
    spectra::AbstractMatrix,
)
    check_compact_grid(m, ν)
    axes(spectra, 1) == Base.OneTo(size(spectra, 1)) || throw(ArgumentError(
        "the input spectra must use one-based row indexing",
    ))
    axes(spectra, 2) == Base.OneTo(size(spectra, 2)) || throw(ArgumentError(
        "the input spectra must use one-based column indexing",
    ))
    axes(output) == (
        Base.OneTo(size(output, 1)),
        Base.OneTo(size(output, 2)),
    ) || throw(ArgumentError(
        "the output spectra must use one-based indexing",
    ))
    size(spectra, 1) == length(ν) || throw(DimensionMismatch(
        "spectra and input grid lengths do not match",
    ))
    size(output) == (length(m.ν_out), size(spectra, 2)) ||
        throw(DimensionMismatch(
            "output dimensions do not match the compact ILS detector grid and spectra",
        ))
    Base.mightalias(output, spectra) && throw(ArgumentError(
        "compact ILS input and output arrays must not alias",
    ))
    @inbounds for spectrum_index in axes(spectra, 2)
        for channel in axes(output, 1)
            first_index = m.first_indices[channel]
            count = m.lengths[channel]
            output[channel, spectrum_index] = dot(
                view(m.weights, 1:count, channel),
                view(
                    spectra,
                    first_index:(first_index + count - 1),
                    spectrum_index,
                ),
            )
        end
    end
    return output
end

function conv_spectra(
    m::CompactVariableKernelInstrument{FT},
    ν::AbstractVector,
    spectrum::AbstractVector{ST},
) where {FT,ST}
    output = Vector{promote_type(FT, ST)}(undef, length(m.ν_out))
    return conv_spectra!(output, m, ν, spectrum)
end

function conv_spectra(
    m::CompactVariableKernelInstrument{FT},
    ν::AbstractVector,
    spectra::AbstractMatrix{ST},
) where {FT,ST}
    output = Matrix{promote_type(FT, ST)}(
        undef,
        length(m.ν_out),
        size(spectra, 2),
    )
    return conv_spectra!(output, m, ν, spectra)
end

spectral_range(ν) = range(first(ν), last(ν); length=length(ν))
