using Distributions
using InstrumentOperator
using Test

@testset "fixed kernel convolution" begin
    ν = collect(-0.01:0.0005:0.01)
    spectrum = ones(length(ν))
    kernel = create_instrument_kernel(Normal(0, 0.0004), -0.002:0.0005:0.002)
    instrument = FixedKernelInstrument(kernel, collect(-0.005:0.001:0.005))
    y = conv_spectra(instrument, ν, spectrum)
    @test length(y) == length(instrument.ν_out)
    @test all(isfinite, y)
    @test all(0.9 .< y .< 1.1)
end

@testset "variable OCO-style kernel preparation" begin
    grid_x = -0.001:0.0005:0.001
    npos = 5
    ils_delta = zeros(3, npos, 1, 1)
    ils_response = zeros(3, npos, 1, 1)
    for i in 1:npos
        ils_delta[:, i, 1, 1] .= [-0.0005, 0.0, 0.0005]
        ils_response[:, i, 1, 1] .= [0.0, 1.0, 0.0]
    end
    kernel = InstrumentOperator.prepare_ils_table(grid_x, ils_response, ils_delta, [1, 1])
    instrument = VariableKernelInstrument(kernel, collect(-0.002:0.001:0.002), collect(1:npos))
    ν = collect(-0.005:0.0005:0.005)
    y = conv_spectra(instrument, ν, ones(length(ν)))
    @test length(y) == npos
    @test all(isfinite, y)
end
