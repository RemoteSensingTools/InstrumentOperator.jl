module InstrumentOperatorCUDAExt

using CUDA
using InstrumentOperator

function InstrumentOperator.gpu_operator(
    instrument::InstrumentOperator.CompactVariableKernelInstrument,
)
    return InstrumentOperator.CompactVariableKernelInstrument(
        CuArray(instrument.weights),
        CuArray(instrument.first_indices),
        CuArray(instrument.lengths),
        CuArray(instrument.ν_in),
        CuArray(instrument.ν_out),
    )
end

function compact_vector_kernel!(
    output,
    weights,
    first_indices,
    lengths,
    spectrum,
    channel_count,
)
    channel = (blockIdx().x - 1) * blockDim().x + threadIdx().x
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
    return
end

function compact_matrix_kernel!(
    output,
    weights,
    first_indices,
    lengths,
    spectra,
    channel_count,
    spectrum_count,
)
    channel = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    spectrum_index = (blockIdx().y - 1) * blockDim().y + threadIdx().y
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
    return
end

function InstrumentOperator.conv_spectra!(
    output::CuVector,
    instrument::InstrumentOperator.CompactVariableKernelInstrument,
    wavelength::CuVector,
    spectrum::CuVector,
)
    InstrumentOperator.check_compact_grid(instrument, wavelength)
    length(spectrum) == length(wavelength) || throw(DimensionMismatch(
        "spectrum and input grid lengths do not match",
    ))
    channel_count = length(instrument.ν_out)
    length(output) == channel_count || throw(DimensionMismatch(
        "output length does not match the compact ILS detector grid",
    ))
    threads = 256
    blocks = cld(channel_count, threads)
    @cuda threads=threads blocks=blocks compact_vector_kernel!(
        output,
        instrument.weights,
        instrument.first_indices,
        instrument.lengths,
        spectrum,
        channel_count,
    )
    return output
end

function InstrumentOperator.conv_spectra!(
    output::CuMatrix,
    instrument::InstrumentOperator.CompactVariableKernelInstrument,
    wavelength::CuVector,
    spectra::CuMatrix,
)
    InstrumentOperator.check_compact_grid(instrument, wavelength)
    size(spectra, 1) == length(wavelength) || throw(DimensionMismatch(
        "spectra and input grid lengths do not match",
    ))
    channel_count = length(instrument.ν_out)
    spectrum_count = size(spectra, 2)
    size(output) == (channel_count, spectrum_count) ||
        throw(DimensionMismatch(
            "output dimensions do not match the detector grid and spectra",
        ))
    threads = (16, 16)
    blocks = (
        cld(channel_count, threads[1]),
        cld(spectrum_count, threads[2]),
    )
    @cuda threads=threads blocks=blocks compact_matrix_kernel!(
        output,
        instrument.weights,
        instrument.first_indices,
        instrument.lengths,
        spectra,
        channel_count,
        spectrum_count,
    )
    return output
end

end
