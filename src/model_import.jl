# used for constructing appropraite dict from YAML object.
function construct_type_map(t::Symbol, constructor::YAML.Constructor,
                            node::YAML.Node)
    mapping = _symbol_dict(YAML.construct_mapping(constructor, node))
    mapping[:tag] = t
    mapping
end

function guess_model_type(data)
    # if the yaml file has the model type key, use what is given there.
    if haskey(data, "model_type")
        return symbol(data["model_type"])
    end

    # othewise we need to do a bit more guesswork
    if haskey(data["symbols"], "shocks")
        if typeof(data["equations"]) == Dict{Any,Any}
            return :dtcscc
        else
            return :dynare
        end
    else
        return :dtmscc
    end
end

function yaml_import(::Type{SymbolicModel}, url; print_code::Bool=false)
    funcs = Dict("!Cartesian" => (c, n) -> construct_type_map(:Cartesian, c, n),
                 "!Normal" => (c, n) -> construct_type_map(:Normal, c, n))

    if match(r"(http|https):.*", url) != nothing
        res = get(url)
        buf = IOBuffer(res.data)
        data = load(buf, funcs)
    else
        data = load_file(url, funcs)
    end

    fname = basename(url)

    model_type = guess_model_type(data)

    SymbolicModel(data, model_type, fname)
end

function yaml_import(url; print_code::Bool=false)

    sm = yaml_import(SymbolicModel, url; print_code=print_code)

    if sm.model_type == :dtcscc
        return DTCSCCModel(sm; print_code=print_code)
    elseif sm.model_type == :dtmscc
        return DTMSCCModel(sm; print_code=print_code)
    else
        throw(Exception)
    end

end
