# Usage examples:
# dict = YAML.load_file("json/oco2.yaml"; dicttype=OrderedDict{String,Any})
# met = Dataset(oco_met_file)
# ocoData = Dataset(oco_file)
# geometry, ils, meteo, measurement = load_L1(dict,ocoData, met)

function load_L1(dict, L1, met)
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

    return geometry, ils, meteo, measurement
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

