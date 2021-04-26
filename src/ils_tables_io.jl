using JSON

function read_ils_table(file::String, jsonFile::String)
    @info "reading ILS file $file" 
    if isfile(file) && isfile(jsonFile)
        jsonDict     = JSON.parsefile(jsonFile);
        ils_json     = jsonDict["ILS"]
        ilsFile      = Dataset(file);
        ils_grid     = getNC_var(ilsFile, ils_json["ils_grid"])
        ils_response = getNC_var(ilsFile, ils_json["ils_response"])
        dispersion   = getNC_var(ilsFile, ils_json["dispersion"])
        @info "ILS table size: " size(ils_response)
        close(ilsFile);
        return ils_grid, ils_response, dispersion
    else
        @error "ILS files don't exist" file isfile(file) jsonFile isfile(jsonFile)
        return nothing, nothing
    end
end

function getNC_var(fin, path)
    try
        loc = split(path, r"/")
        # println(loc)
        if length(loc) == 1
            return fin[path].var[:]
        elseif length(loc) > 1
            gr = []
            for i in 1:length(loc) - 1
                if i == 1
                    gr = fin.group[loc[i]]
                else
                    gr = gr.group[loc[i]]
                end
            end
            return gr[loc[end]].var[:]
        end
    catch e
        @error e
        println("Error in getNC_var ", path)
        return nothing
    end
end