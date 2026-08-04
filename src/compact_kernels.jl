@kernel function compact_vector_kernel!(
    output,
    weights,
    first_indices,
    lengths,
    spectrum,
    channel_count,
)
    channel = @index(Global, Linear)
    if channel <= channel_count
        first_index = first_indices[channel]
        count = lengths[channel]
        value = zero(eltype(output))
        @inbounds for local_index in 1:count
            value += weights[local_index, channel] *
                spectrum[first_index + local_index - 1]
        end
        @inbounds output[channel] = value
    end
end

@kernel function compact_matrix_kernel!(
    output,
    weights,
    first_indices,
    lengths,
    spectra,
    channel_count,
    spectrum_count,
)
    index = @index(Global, Cartesian)
    channel = index[1]
    spectrum_index = index[2]
    if channel <= channel_count && spectrum_index <= spectrum_count
        first_index = first_indices[channel]
        count = lengths[channel]
        value = zero(eltype(output))
        @inbounds for local_index in 1:count
            value += weights[local_index, channel] *
                spectra[first_index + local_index - 1, spectrum_index]
        end
        @inbounds output[channel, spectrum_index] = value
    end
end
