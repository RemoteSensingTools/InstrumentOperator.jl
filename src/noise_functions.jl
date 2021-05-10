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
    @unpack t_int, Δλ, Qₑ, η = ins
    # Common unit:
    #luminosity = uconvert.(u"J/s/μm^2/sr/μm",radiance)
    energy = radiance * etendue(ins) * t_int * Δλ * Qₑ * η
    uconvert.(u"nm^-1", energy, Spectral()) .* λ
end

function noise_at_fpa(ins, photons)
    @unpack t_int, dark_current, σ_read = ins
    total_electrons = photons .+ dark_current .* t_int
    # Sum variances; Can add more noise terms here later if needed!
    total_variance = total_electrons .+ σ_read^2 
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

