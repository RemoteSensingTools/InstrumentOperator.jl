module InstrumentOperatorCUDAExt

using CUDA
using InstrumentOperator
using KernelAbstractions

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
    backend = KernelAbstractions.get_backend(output)
    kernel! = InstrumentOperator.compact_vector_kernel!(backend, 256)
    kernel!(
        output,
        instrument.weights,
        instrument.first_indices,
        instrument.lengths,
        spectrum,
        channel_count,
        ndrange=channel_count,
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
    backend = KernelAbstractions.get_backend(output)
    kernel! = InstrumentOperator.compact_matrix_kernel!(backend, (16, 16))
    kernel!(
        output,
        instrument.weights,
        instrument.first_indices,
        instrument.lengths,
        spectra,
        channel_count,
        spectrum_count,
        ndrange=(channel_count, spectrum_count),
    )
    return output
end

end
