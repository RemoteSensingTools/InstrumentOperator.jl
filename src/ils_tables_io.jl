using JSON

function read_ils_table(file::String, jsonFile::String)
    @info "reading ILS file $file" 
    isfile(file) || throw(ArgumentError("ILS file does not exist: $file"))
    isfile(jsonFile) || throw(ArgumentError(
        "ILS metadata file does not exist: $jsonFile",
    ))
    jsonDict = JSON.parsefile(jsonFile)::Dict{String,Any}
    ils_json = jsonDict["ILS"]::Dict{String,Any}
    ilsFile = Dataset(file)
    try
        ils_grid = getNC_var(ilsFile, ils_json["ils_grid"]::String)
        ils_response = getNC_var(ilsFile, ils_json["ils_response"]::String)
        dispersion = getNC_var(ilsFile, ils_json["dispersion"]::String)
        @info "ILS table size: " size(ils_response)
        return ils_grid, ils_response, dispersion
    finally
        close(ilsFile)
    end
end

function getNC_var(fin, path)
    loc = split(path, r"/")
    if length(loc) == 1
        return fin[path].var
    end
    gr = fin.group[first(loc)]
    for i in 2:(length(loc) - 1)
        gr = gr.group[loc[i]]
    end
    return gr[last(loc)].var
end
