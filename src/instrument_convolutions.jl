
"Convolves and resamples the input spectrum with a fixed kernel"
function conv_spectra(m::FixedKernelInstrument, ν, spectrum)
    s = imfilter(spectrum, m.kernel)
    interp_cubic = CubicSplineInterpolation(ν, s)
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
    te = LinearInterpolation(m.ν_out, FT.(m.ind_out); extrapolation_bc=Interpolations.Flat())
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
    fin = LinearInterpolation(ν[ind], spec_out; extrapolation_bc=Interpolations.Flat())
    return fin(m.ν_out)
end;


