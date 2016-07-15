# used for constructing appropraite dict from YAML object.
function construct_type_map(t::Symbol, constructor::YAML.Constructor,
                            node::YAML.Node)
    mapping = _symbol_dict(YAML.construct_mapping(constructor, node))
    mapping[:tag] = t
    mapping
end

function guess_model_type(data)
    # if the yaml file has the model type key, use what is given there.
    if haskey(data, :model_type)
        return Symbol(data[:model_type])
    end

    # othewise we need to do a bit more guesswork
    if haskey(data[:symbols], :shocks)
        if typeof(data[:equations]) == Dict{Any,Any}
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
        data = _symbol_dict(load(buf, funcs))
    else
        data = _symbol_dict(load_file(url, funcs))
    end

    fname = basename(url)

    model_type = guess_model_type(data)

    SymbolicModel(data, model_type, fname)
end

function yaml_import(url; print_code::Bool=false)
    sm = yaml_import(SymbolicModel, url; print_code=print_code)
    NumericModel(sm; print_code=print_code)
end

function dynare_import(::Type{SymbolicModel}, url; print_code::Bool=false)
    funcs = Dict("!Cartesian" => (c, n) -> construct_type_map(:Cartesian, c, n),
                 "!Normal" => (c, n) -> construct_type_map(:Normal, c, n))

    if match(r"(http|https):.*", url) != nothing
        res = get(url)
        buf = IOBuffer(res.data)
        data = dynare_parser(readlines(buf), basename(url))
    else
        data = load_modfile(url)
    end
end

function dynare_import(url; print_code::Bool=false)
    sm = dynare_import(SymbolicModel, url; print_code=print_code)
    m = NumericModel(sm; print_code=print_code)
    # also need to compile first and second derivative for dynare model
    make_method(Der{1}, m.factories[:dynare])
    make_method(Der{2}, m.factories[:dynare], mutating=false)
    m
end

function _extract_calib_block(text, regex)
    out = OrderedDict{Symbol,Union{Expr,Symbol,Number}}()
    tmp  = match(regex, text)
    if tmp != nothing  # if we matched
        tmp = split(tmp[1], ";")
        filter!(x -> ismatch(r"\S", x), tmp) # Remove lines if no chars present
        for ln in tmp
            pair = match(r"(.+)=(.+)", ln)
            out[Symbol(strip(pair[1]))] = parse(strip(pair[2]))
        end
    end
    out
end

function _extract_variable_group(text, regex)
    tmp = match(regex, text)
    out = split(tmp[1], " ")
    filter!(x -> ismatch(r"\S", x), out) # Remove lines if no chars present
    map(Symbol, out)
end

function load_modfile(modfile_name::String)
    lines = open(readlines, modfile_name)
    dynare_parser(lines, modfile_name)
end

function dynare_parser(lines::Vector, modfile_name="nofile")
    # Remove lines beginning with comments (%, //, /*, \*)
    filter!(x -> !ismatch(r"^\s*%", x), lines)
    filter!(x -> !ismatch(r"^\s*\/\/", x), lines)
    filter!(x -> !ismatch(r"^\s*\/\*", x), lines)

    # Remove end-of-line comments
    for ln = 1:length(lines)
        lines[ln] = split(lines[ln], r"%.*$")[1]    # %  comments
        lines[ln] = split(lines[ln], r"\/\/.*$")[1]  # \\ comments
        lines[ln] = split(lines[ln], r"/\*.*$")[1]   # \* comments
    end

    # Smoosh text back together
    text = join(lines)

    # Remove new lines, spaces, carriage returns,
    text = replace(text, "\t", " ")
    text = replace(text, "\r", " ")
    text = replace(text, "\n", " ")

    # Get variable names
    variables = _extract_variable_group(text, r"var\s(.*?);")
    shocks = _extract_variable_group(text, r"varexo\s(.*?);")
    parameters = _extract_variable_group(text, r"parameters\s(.*?);")

    # get equations, fill a list
    tmp  = match(r"model;\s(.*?)end;", text)
    equations = split(tmp[1], ";")
    filter!(x -> ismatch(r"\S", x), equations) # Remove lines if no chars present
    equations = map(parse, equations)

    # Get calibration values, fill some dicts
    param_values = _extract_calib_block(text, r"parameters\s(?:.*?);(.*)model")
    initval = _extract_calib_block(text, r"initval;(.*?)end;")
    endval = _extract_calib_block(text, r"endval;(.*?)end;")

    # Get the calibrated shock values, fill matrix
    shockvaldict = Dict()
    shock_matrix = Array(Any, length(shocks), length(shocks))
    fill!(shock_matrix, 0)
    tmp = match(r"shocks;(.*?)(.*?)end;", text)
    if tmp != nothing
        if contains(tmp[2], "stderr")
            tmpkey = matchall(r"var(.*?);", tmp[2])
            tmpentry = matchall(r"stderr(.*?);", tmp[2])
            for ln = 1:length(tmpkey)
                key = strip(match(r"var\s(.*)", tmpkey[ln])[1])
                key = match(r"(.*);", key)[1]
                entry = strip(match(r"stderr\s(.*)", tmpentry[ln])[1])
                entry = match(r"(.*);", entry)[1]
                shockvaldict[key] =  entry
            end
        else
            tmp = split(tmp[2], ";")
            filter!(x -> ismatch(r"\S", x), tmp) # Remove lines if no chars present
            for ln = 1:length(tmp)
                key = match(r"var\s(.*)\s=", tmp[ln])[1]
                entry = match(r"=\s(.*)", tmp[ln])[1]
                shockvaldict[key] =  entry
            end
        end
        # Fill matrix in order of shocks in "shocks" dictionary
        for ln = 1:length(shocks)
            if haskey(shockvaldict, shocks[ln])
                shock_matrix[ln, ln] = parse(shockvaldict[shocks[ln]])
            else
                shock_matrix[ln, ln] = 0.0
            end
        end
    end

    symbols = OrderedDict{Symbol,Vector{Symbol}}(:variables => variables,
                                                :shocks => shocks,
                                                :parameters => parameters)
    eqs = OrderedDict{Symbol,Vector{Expr}}(:dynare => equations)
    calib = merge(param_values, initval)
    distribution = Dict(:tag => :Normal, :sigma => shock_matrix)
    options = Dict{Symbol,Any}(:distribution => distribution, :endval => endval)
    defs = Dict()
    recipe = RECIPES[:dynare]

    # go through parameters, shocks, and variables and make sure they are all
    # present in calib. If they aren't, follow the dynare convention and add an
    # entry in calib that sets them to 0
    for (_, syms) in symbols
        for s in syms
            get!(calib, s, 0.0)
        end
    end

    _name = string(split(basename(modfile_name), ".mod")[1])
    _id = gensym(Symbol(_name))

    SymbolicModel{_id,:dynare}(recipe, symbols, eqs, calib, options, defs,
                               _name, modfile_name)
end
