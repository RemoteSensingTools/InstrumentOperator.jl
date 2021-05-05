"Compute etendue of the instrument"
function etendue(instrument::GratingNoiseModel)
    @unpack detector_size, Fnumber = instrument
    "Compute detector area"
    A = detector_size^2;
    "Compute solid angle Ω"
    Ω = uconvert(u"sr",π/(4Fnumber^2))
    "Compute etendue (same unit as area A)"
    A * Ω
end

"Compute photons converted to electrons at the FPA"
function photons_at_fpa(ins, λ, radiance)
    @unpack integration_time, SSI, FPA_qe, grating_efficiency, effTransmission = ins
    # Common unit:
    luminosity = uconvert.(u"J/s/μm^2/sr/μm",radiance)
    total_efficiency = FPA_qe * grating_efficiency * effTransmission
    energy = luminosity * etendue(ins) * integration_time * SSI * total_efficiency
    uconvert.(u"nm^-1", energy, Spectral()) .* λ
end

function noise_at_fpa(ins, photons)
    @unpack integration_time, dark_current, readNoise = ins
    total_electrons = photons .+ dark_current .* integration_time
    # shot noise, variance = electrons
    shot_noise_variance = total_electrons
    # Sum variances; Can add more noise terms here later if needed!
    total_variance = shot_noise_variance .+ readNoise^2
    # Return 1-σ noise level
    sqrt.(total_variance)
end

"Compute noise equivalent spectral radiance"
function noise_equivalent_radiance(ins, λ, radiance)
    pfpa  = photons_at_fpa(ins, λ, radiance)
    noise = noise_at_fpa(ins,pfpa)
    SNR = pfpa./noise
    return radiance./SNR
end

