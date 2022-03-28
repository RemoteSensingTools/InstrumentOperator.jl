# Usage examples:
# dict = YAML.load_file("json/oco2.yaml"; dicttype=OrderedDict{String,Any})
# met = Dataset(oco_met_file)
# ocoData = Dataset(oco_file)
# geometry, ils, meteo, measurement = load_L1(dict,ocoData, met)

function load_L1(dict, L1::NCDataset, met::NCDataset)
    # Dicts with mostly NCVars but also some basic stuff:
    geometry    = Dict{String, Any}()
    ils         = Dict{String, Any}()
    meteo       = Dict{String, Any}()
    measurement = Dict{String, Any}()
    
    # Reading in Geometries
    loadEntries!(L1, dict, geometry, "Geometry")

    # Reading in Meterorology
    loadEntries!(met, dict, meteo, "Meteorology")
    
    # Reading in Meausurements
    loadEntries!(L1, dict, measurement, "Measurement")
    
    # Reading in Meausurements
    loadEntries!(L1, dict, ils, "ILS")

    return L1_OCO(geometry, ils, meteo, measurement)
end

function load_L1(dict::String, L1::String, met::String)
    ocoData = Dataset(L1);    
    metData = Dataset(met);
    # Load dictionary:
    dictOCO2 = YAML.load_file(dict; dicttype=OrderedDict{String,Any});

    # Load L1 file (could just use filenames here as well later on)
    return load_L1(dictOCO2,ocoData, metData)
end

function loadEntries!(NC, dict, out_dict, name)
    if haskey(dict,name)
        for (key, value) in dict[name]
            if typeof(value) == String
                out_dict[key] = getNC_var(NC, value)
            else
                out_dict[key] = value
            end
        end
    end
end

function doppler_factor(vᵣ::FT) where FT
    c = FT(299792458)
    sqrt((c-vᵣ)/(c+ vᵣ))
end


function getMeasurement(oco::L1_OCO, bands::Tuple, indices::Tuple, GeoInd; kernel_range = 0.45e-3,kernel_step = 0.001*1e-3 )
    @assert length(indices) == length(bands) "Length of bands and indices has to be identical"
    n = length(indices)
    # First band:
    band = oco.measurement["bands"][bands[1]]
    # Extended Dimensions (irrespective of sounding)
    extended_dims = [GeoInd[1],bands[1]];
    ind  = indices[1]
    bandIndices = (ind .- ind[1] .+ 1,)
    rad = oco.measurement[band][ind,GeoInd...]
    FT = typeof(rad[1])
    dispPoly = Polynomial(view(oco.ils["dispersion"], :, extended_dims...))
    f_doppler = doppler_factor(oco.geometry["v_rel"][GeoInd[2]]);
    # Apply doppler shift (all depends on definitions)
    # @info("Doppler Shift factor = $(f_doppler)")
    ν = FT.(dispPoly.(indices[1])) 
    # First ILS
    # First hard-coded:
    #@show FT, typeof(ν)
    
    # Hard coded for now, needs to be changed later:
    grid_x = FT(-kernel_range):FT(kernel_step):FT(kernel_range)
    ils_pixel   = prepare_ils_table(grid_x, oco.ils["ils_response"][:], oco.ils["ils_grid"][:],extended_dims)
    oco2_kernels = (VariableKernelInstrument(ils_pixel, ν, collect(ind .-1)),)
    # Concatenate rest (if applicable)
    for i=2:n
        band = oco.measurement["bands"][bands[i]]
        ind  = indices[i]
        bandIndices = (bandIndices..., ind .-ind[1] .+ 1 .+ bandIndices[end][end])
        rad = [rad; oco.measurement[band][ind,GeoInd...]]
        extended_dims = [GeoInd[1],bands[i]];
        dispPoly = Polynomial(view(oco.ils["dispersion"], :, extended_dims...))
        ν = [ν; FT.(dispPoly.(indices[i]))]
        # ILS kernels:
        # grid_x = FT(-0.35e-3):FT(0.001*1e-3):FT(0.35e-3)
        ils_pixel   = prepare_ils_table(grid_x, oco.ils["ils_response"][:], oco.ils["ils_grid"][:],extended_dims)
        #@show ind
        oco2_kernels = (oco2_kernels..., VariableKernelInstrument(ils_pixel, FT.(dispPoly.(indices[i])), collect(ind .-1)))
        #@show VariableKernelInstrument(ils_pixel, FT.(dispPoly.(indices[i])), collect(ind .-1)).ν_out
    end

    #### Meteo stuff  ###
    p_half = reverse(oco.meteo["ak"] + oco.meteo["bk"] * oco.meteo["p_surf"][GeoInd...]/100)
    T   =    oco.meteo["T"][:,GeoInd...] 
    q   =    oco.meteo["q"][:,GeoInd...]
    lat = oco.geometry["lat"][GeoInd...]
    lon = oco.geometry["lon"][GeoInd...]
    sza = oco.geometry["sza"][GeoInd...]
    ϕ   = oco.geometry["ϕ"][GeoInd...]
    vza = oco.geometry["vza"][GeoInd...]
    azi = oco.geometry["azi"][GeoInd...]

    # From OCO-2 ATBD Page 53 https://docserver.gesdisc.eosdis.nasa.gov/public/project/OCO/OCO_L1B_ATBD.pdf
    StokesCoef = [FT(0.5),FT(cosd(2ϕ)/2), FT(sind(2ϕ)/2), FT(0)] 

    return MeasurementOCO(
        ν,
        bandIndices,
        rad,
        oco.measurement["radiance_o2"].attrib["Units"],
        "μm",
        lat,
        lon,
        vza,
        sza,
        azi,
        FT((1 + 0.0053024*sind(lat)^2 - 0.0000058*sin(2lat)^2) * 9.780318), # WELMEC formula
        ϕ,
        StokesCoef,
        FT.(p_half),
        FT.((p_half[1:end-1] + p_half[2:end])/2),
        T,
        q,
        oco2_kernels,
        f_doppler
    )
end

