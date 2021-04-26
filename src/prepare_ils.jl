
"""
$(FUNCTIONNAME)(grid_x::AbstractRange, ils_in::Array{FT}, ils_Δ::Array{FT}, extended_dims::Array{Int}=[]) where {FT <: AbstractFloat}

Pre-compute the ILS table input as function of spectral distance from center converted to the modeling grid 
Input: grid_x, ils_in, ils_Δ, extended_dims
Output: Offset-Array with tabulated responses interpolated to grid_x
"""
function prepare_ils_table(grid_x::AbstractRange, ils_in::Array{FT}, ils_Δ::Array{FT}, extended_dims::Array{Int}=[]) where {FT <: AbstractFloat}
    @assert minimum(abs.(grid_x)) < eps(FT) "grid_x must include 0 (center pixel)"
    ind_0 = argmin(abs.(grid_x))
    axis_pixel = (-ind_0 + 1):(grid_x.len - ind_0)
    # Dimension of ILS table (at least 1D, first dimension needs to be across wavenumber/wavelength)
    dims  = size(ils_in);
    # number of spectral positions of ILS table
    n_x = dims[1]
    # Number of ILS per detector position (can be 1 if ILS is constant across detector grid)
    n_pos = dims[2]
    ils    = view(ils_in, :, :,  extended_dims...);
    ils_Δ_ = view(ils_Δ,  :, :,  extended_dims...);
    ils_pixel = zeros(FT, grid_x.len, n_pos);

    for i = 1:n_pos
        ind = findall(minimum(ils_Δ_[:,i]) .< grid_x .< maximum(ils_Δ_[:,i]));
        interp = Interpolations.LinearInterpolation(ils_Δ_[:,i], ils[:,i])
        ils_pixel[ind,i] = interp.(grid_x[ind]);
    end
    # normalize here
    return OffsetArray(ils_pixel ./ sum(ils_pixel, dims=1), axis_pixel, 1:n_pos)
end
# This can still derive min/max range automatically, need to double check

function create_instrument_kernel(FTS::FTSInstrument, grid_x::AbstractRange, ν̅)
    @unpack MOPD, β, FOV = FTS
    # Need to make sure this is 0!
    i₀ = argmin(abs.(grid_x))
    axis_pixel = (-i₀ + 1):(grid_x.len - i₀)
    @info "Center of OffsetArrays is " minimum(abs.(grid_x))
    #Δσ = 1/(2FTS.OPD)
    
    # Create sinc FTS kernel
    XX = 2MOPD * grid_x
    sinc_kernel = OffsetArray(2MOPD * sinc.(XX),axis_pixel)
    
    # Now to FOV box:
    Θ = ν̅/2 * FOV^2
    β = β * Θ
    @show Θ, β
    # Create empty box_kernel:
    box = collect(0*grid_x);
    for i in eachindex(grid_x)
        if grid_x[i] ≤ 0
            if Θ > β && grid_x[i] ≥ -(Θ - β)
                box[i] = 1
            elseif grid_x[i] ≥ -(Θ + β)
                box[i] = acos((grid_x[i]^2 + β^2 - Θ^2)/(2grid_x[i] * β)) / π
            end
        end
    end
    box_kernel = OffsetArray(box,axis_pixel)
    box_kernel /= sum(box_kernel)
    sinc_kernel /= sum(sinc_kernel)
    return imfilter(box_kernel,sinc_kernel)
end

"Create kernel from a single continous distribution function"
function create_instrument_kernel(di::ContinuousUnivariateDistribution, grid_x::AbstractRange)
    ils = pdf.(di,grid_x)
    i₀ = argmin(abs.(grid_x))
    axis_pixel = (-i₀ + 1):(grid_x.len - i₀)
    ils /= sum(ils)
    return OffsetArray(ils,axis_pixel)
end

"Create kernel from multiple single continous distribution function (convolves those)"
function create_instrument_kernel(di::Array{ContinuousUnivariateDistribution}, grid_x::AbstractRange)
    ils = pdf.(di[1],grid_x)
    for i=2:length(di)
        ils = imfilter(ils, pdf.(di[i],grid_x))
    end
    i₀ = argmin(abs.(grid_x))
    axis_pixel = (-i₀ + 1):(grid_x.len - i₀)
    ils /= sum(ils)
    return OffsetArray(ils,axis_pixel)
end