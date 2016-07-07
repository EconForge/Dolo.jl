Normal(;sigma=zeros(0, 0)) = Normal(sigma)
Cartesian(;a=[], b=[], orders=[]) = Cartesian(a, b, orders)

function construct_type_map(t::Symbol, constructor::YAML.Constructor, node::YAML.Node)
    mapping = _symbol_dict(YAML.construct_mapping(constructor, node))
    mapping[:kind] = t
    mapping
end

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

function yaml_import(url; print_code::Bool=false)

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

    sym_model = SymbolicModel(data, model_type, fname)

    if model_type == :dtcscc
        return DTCSCCModel(sym_model; print_code=print_code)
    elseif model_type == :dtmscc
        return DTMSCCModel(sym_model; print_code=print_code)
    else
        throw(Exception)
    end

end
