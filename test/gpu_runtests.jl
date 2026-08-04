using CUDA
using InstrumentOperator
using Test

CUDA.functional() || error("a functional CUDA device is required")
CUDA.allowscalar(false)

wavelength = Float32.(collect(757.68:0.001:758.32))
detector = Float32.(collect(757.75:0.02:758.25))
delta = Float32.(collect(-0.06:0.002:0.06))
response = repeat(
    reshape(exp.(-0.5f0 .* (delta ./ 0.025f0).^2), :, 1),
    1,
    length(detector),
)
offset = repeat(reshape(delta, :, 1), 1, length(detector))
cpu_operator = prepare_compact_ils(
    wavelength,
    detector,
    response,
    offset,
)
gpu = gpu_operator(cpu_operator)

spectra = Float32.(
    sin.(wavelength .* 20) .+
    reshape(collect(range(0, 1; length=128)), 1, :)
)
cpu = conv_spectra(cpu_operator, wavelength, spectra)
device_spectra = CuArray(spectra)
CUDA.@sync device = conv_spectra(gpu, gpu.ν_in, device_spectra)

@test Array(device) ≈ cpu rtol=2f-6 atol=2f-6
@test maximum(abs.(Array(sum(gpu.weights; dims=1)) .- 1f0)) < 2f-6
