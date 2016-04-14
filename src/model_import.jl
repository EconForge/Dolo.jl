
function guess_model_type(data)
    if ("shocks" in keys(data["symbols"]))
        if typeof(data["equations"]) == Dict{Any,Any}
            return :dtcscc
        else
            return :dynare
        end
    else
        return :dtmscc
    end
end

"""
Imports the model from a yaml file specified by the `url` input
parameter, and returns the corresponding `Dolo` model object.
"""
function yaml_import(url)

    if match(r"(http|https):.*", url) != nothing
        res = get(url)
        buf = IOBuffer(res.data)
        data = load(buf)
    else
        data = load_file(url)
    end

    fname = basename(url)

    model_type = guess_model_type(data)

    sym_model = SymbolicModel(data, model_type, fname)

    if model_type == :dtcscc
        return DTCSCCModel(sym_model)
    elseif model_type == :dtmscc
        return DTMSCCModel(sym_model)
    else
        throw(Exception)
    end

end
