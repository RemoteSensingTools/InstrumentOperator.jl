using Aqua
using Distributions
using InstrumentOperator
using JET
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

@testset "quality tooling" begin
    Aqua.test_all(InstrumentOperator)
    if get(ENV, "INSTRUMENT_OPERATOR_RUN_JET", "false") == "true"
        JET.test_package(
            InstrumentOperator;
            target_modules=(InstrumentOperator,),
        )
    end
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

@testset "compact exact variable ILS" begin
    FT = Float32
    ν_in = FT.(collect(-0.01:0.0001:0.01))
    ν_out = FT[-0.002, 0.0, 0.002]
    delta = repeat(
        reshape(FT.(collect(-0.0006:0.0002:0.0006)), :, 1),
        1,
        length(ν_out),
    )
    response = similar(delta)
    for channel in axes(response, 2)
        width = FT(0.00022 + 0.00001 * channel)
        response[:, channel] .= exp.(-FT(0.5) .* (
            delta[:, channel] ./ width
        ).^2)
    end

    instrument = prepare_compact_ils(ν_in, ν_out, response, delta)
    @test instrument isa CompactVariableKernelInstrument{Float32}
    @test all(instrument.lengths .> 2)
    for channel in eachindex(instrument.lengths)
        @test sum(view(
            instrument.weights,
            1:instrument.lengths[channel],
            channel,
        )) ≈ 1f0 atol=2f-7
    end

    flat = ones(FT, length(ν_in))
    flat_output = conv_spectra(instrument, ν_in, flat)
    @test flat_output ≈ ones(FT, length(ν_out)) atol=2f-7

    linear = FT.(1 .+ 2 .* ν_in)
    linear_output = conv_spectra(instrument, ν_in, linear)
    @test linear_output ≈ FT.(1 .+ 2 .* ν_out) atol=3f-5

    spectra = hcat(flat, linear, linear .^ 2)
    matrix_output = conv_spectra(instrument, ν_in, spectra)
    @test size(matrix_output) == (length(ν_out), 3)
    for column in axes(spectra, 2)
        @test matrix_output[:, column] == conv_spectra(
            instrument,
            ν_in,
            spectra[:, column],
        )
    end

    output = similar(flat_output)
    conv_spectra!(output, instrument, ν_in, flat)
    @test output == flat_output
    @test (@allocated conv_spectra!(output, instrument, ν_in, flat)) == 0

    double_output = conv_spectra(instrument, ν_in, Float64.(linear))
    @test eltype(double_output) == Float64
    @test double_output ≈ Float64.(linear_output) atol=2e-7

    shifted_grid = copy(ν_in)
    shifted_grid[10] += FT(1e-6)
    @test_throws ArgumentError conv_spectra(
        instrument,
        shifted_grid,
        flat,
    )
    @test_throws ArgumentError prepare_compact_ils(
        ν_in,
        ν_out,
        response,
        delta;
        half_width=FT(-0.001),
    )
end
