using CUDA
using InstrumentOperator
using Printf
using Statistics

CUDA.functional() || error("a functional CUDA device is required")
CUDA.allowscalar(false)

const FT = Float32
const WAVELENGTH = FT.(collect(757.65:0.001:759.20))
const DETECTOR = FT.(collect(757.72:0.014:759.13))
const DELTA = FT.(collect(-0.07:0.002:0.07))
const RESPONSE = repeat(
    reshape(exp.(-FT(0.5) .* (DELTA ./ FT(0.025)).^2), :, 1),
    1,
    length(DETECTOR),
)
const OFFSET = repeat(reshape(DELTA, :, 1), 1, length(DETECTOR))

function median_time(callable; samples=7)
    callable()
    values = [@elapsed callable() for _ in 1:samples]
    return median(values)
end

function main()
    cpu_operator = prepare_compact_ils(
        WAVELENGTH,
        DETECTOR,
        RESPONSE,
        OFFSET,
    )
    device_operator = gpu_operator(cpu_operator)
    @printf(
        "grid=%d channels=%d maximum ILS points=%d\n",
        length(WAVELENGTH),
        length(DETECTOR),
        maximum(cpu_operator.lengths),
    )
    @printf(
        "%8s %12s %12s %12s %14s\n",
        "spectra",
        "CPU ms",
        "GPU ms",
        "speedup",
        "GPU spectra/s",
    )
    for spectrum_count in (1, 8, 128, 1024, 8192, 32768)
        phase_values = spectrum_count == 1 ?
            zeros(FT, 1) :
            FT.(range(0, 2; length=spectrum_count))
        phase = reshape(phase_values, 1, :)
        spectra = @. sin(FT(20) * WAVELENGTH) + phase
        device_spectra = CuArray(spectra)
        cpu_time = spectrum_count <= 8192 ? median_time(
            () -> conv_spectra(cpu_operator, WAVELENGTH, spectra);
            samples=3,
        ) : NaN
        gpu_call = () -> begin
            CUDA.@sync begin
                conv_spectra(
                    device_operator,
                    device_operator.ν_in,
                    device_spectra,
                )
            end
        end
        gpu_time = median_time(gpu_call)
        @printf(
            "%8d %12.3f %12.3f %12.1f %14.0f\n",
            spectrum_count,
            1e3 * cpu_time,
            1e3 * gpu_time,
            cpu_time / gpu_time,
            spectrum_count / gpu_time,
        )
    end
end

main()
