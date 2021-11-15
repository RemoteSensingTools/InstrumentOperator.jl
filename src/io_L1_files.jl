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

function load_L1(dict::OrderedDict, L1::String, met::String)
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

function getMeasurement(oco::L1_OCO, bands::Tuple, indices::Tuple, GeoInd)
    @assert length(indices) == length(bands) "Length of bands and indices has to be identical"
    n = length(indices)
    # First band:
    band = oco.measurement["bands"][bands[1]]
    # Extended Dimensions (irrespective of sounding)
    extended_dims = [GeoInd[1],bands[1]];
    ind  = indices[1]
    bandIndices = (ind .- ind[1] .+ 1,)
    rad = oco.measurement[band][ind,GeoInd...]
    dispPoly = Polynomial(view(oco.ils["dispersion"], :, extended_dims...))
    ν = (dispPoly.(indices[1]))

    # Concatenate rest (if applicable)
    for i=2:n
        band = oco.measurement["bands"][bands[i]]
        ind  = indices[i]
        bandIndices = (bandIndices..., ind .-ind[1] .+ 1 .+ bandIndices[end][end])
        rad = [rad; oco.measurement[band][ind,GeoInd...]]
        extended_dims = [GeoInd[1],bands[i]];
        dispPoly = Polynomial(view(oco.ils["dispersion"], :, extended_dims...))
        ν = [ν; (dispPoly.(indices[i]))]
    end

    #### Meteo stuff  ###
    p_half = reverse(oco.meteo["ak"] + oco.meteo["bk"] * oco.meteo["p_surf"][GeoInd...])/100
    T   =    oco.meteo["T"][:,GeoInd...] 
    q   =    oco.meteo["q"][:,GeoInd...]
    lat = oco.geometry["lat"][GeoInd...]
    lon = oco.geometry["lon"][GeoInd...]
    sza = oco.geometry["sza"][GeoInd...]
    ϕ   = oco.geometry["ϕ"][GeoInd...]
    vza = oco.geometry["vza"][GeoInd...]
    azi = oco.geometry["azi"][GeoInd...]

    return ν, rad, bandIndices, p_half, T, q
end

