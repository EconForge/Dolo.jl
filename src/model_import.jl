# used for constructing appropraite dict from YAML object.
function construct_type_map(t::Symbol, constructor::YAML.Constructor,
                            node::YAML.Node)
    mapping = _symbol_dict(YAML.construct_mapping(constructor, node))
    mapping[:tag] = t
    mapping
end

const yaml_types = let
    pairs = [("!Cartesian", :Cartesian),
             ("!Normal", :Normal),
             ("!MarkovChain", :MarkovChain)]
    Dict{String,Function}([t => (c, n) -> construct_type_map(s, c, n)
                           for (t, s) in pairs])
end

function guess_model_type(data)
    # if the yaml file has the model type key, use what is given there.
    if haskey(data, :model_type)
        return Symbol(data[:model_type])
    end

    # othewise we need to do a bit more guesswork
    haskey(data[:symbols], :shocks) ? :dtcscc : :dtmscc
end

function yaml_import(::Type{SymbolicModel}, url; print_code::Bool=false)
    if match(r"(http|https):.*", url) != nothing
        res = get(url)
        buf = IOBuffer(res.data)
        data = _symbol_dict(load(buf, yaml_types))
    else
        data = _symbol_dict(load_file(url, yaml_types))
    end

    fname = basename(url)

    model_type = guess_model_type(data)

    SymbolicModel(data, model_type, fname)
end

function yaml_import(url; print_code::Bool=false)
    sm = yaml_import(SymbolicModel, url; print_code=print_code)
    NumericModel(sm; print_code=print_code)
end
