# used for constructing appropraite dict from YAML object.
function construct_type_map(t::Symbol, constructor::YAML.Constructor,
                            node::YAML.Node)
    mapping = _symbol_dict(YAML.construct_mapping(constructor, node))
    mapping[:tag] = t
    mapping
end

const yaml_types = let
    pairs = [("!Cartesian", :Cartesian),
            ("!Smolyak", :Smolyak),
             ("!Normal", :Normal),
             ("!MarkovChain", :MarkovChain),
             ("!VAR1", :VAR1)]
    Dict{AbstractString,Function}([(t, (c, n) -> construct_type_map(s, c, n))
                           for (t, s) in pairs])
end



type SModel{ID} <: ASModel{ID}
    data::Dict{Any,Any}
end


function SModel(url::AbstractString)
    if match(r"(http|https):.*", url) != nothing
        res = get(url)
        buf = Dolo.IOBuffer(res.data)
        data = Dolo._symbol_dict(Dolo.load(buf, yaml_types))
    else
        data = Dolo._symbol_dict(Dolo.load_file(url, yaml_types))
    end
    id = gensym()
    return SModel{id}(data)
end

function get_symbols(model::ASModel)
    syms = model.data[:symbols]
    symbols = OrderedDict{Symbol,Vector{Symbol}}()
    for k in keys(syms)
        symbols[Symbol(k)] = [Symbol(e) for e in syms[k]]
    end
    return symbols
end

function get_variables(model::ASModel)
    symbols = get_symbols(model)
    vars = cat(1, values(symbols)...)
    dynvars = setdiff(vars, symbols[:parameters] )
end

function get_name(model::ASModel)
    get(model.data, :name, "modeldoesnotwork")
end

function get_equations(model::ASModel)

    eqs = model.data[:equations]
    recipe = Dolo.RECIPES[:dtcc]
    _symbols = get_symbols(model)

    # prep equations: parse to Expr
    _eqs = OrderedDict{Symbol,Vector{Expr}}()
    for k in keys(recipe[:specs])
        if k in keys(eqs)
            # we handle these separately
            (k in [:arbitrage,]) && continue
            these_eq = get(eqs, (k), [])
            # verify that we have at least 1 equation if section is required
            if !get(recipe[:specs][k], :optional, false)
                length(these_eq) == 0 && error("equation section $k required")
            end
            # finally pass in the expressions
            _eqs[k] = Expr[Dolo._to_expr(eq) for eq in these_eq]
        end
    end
    # handle the arbitrage, arbitrage_exp, controls_lb, and controls_ub
    if haskey(recipe[:specs], :arbitrage) && (:arbitrage in keys(eqs))
        c_lb, c_ub, arb = Dolo._handle_arbitrage(eqs[:arbitrage],
                                            _symbols[:controls])
        _eqs[:arbitrage] = arb
        _eqs[:controls_lb] = c_lb
        _eqs[:controls_ub] = c_ub
    end

    defs = get_definitions(model)
    dynvars = cat(1, get_variables(model), keys(defs)...)

    for eqtype in keys(_eqs)
        _eqs[eqtype] = [sanitize(eq,dynvars) for eq in _eqs[eqtype]]
    end

    return _eqs
end


function get_infos(model::ASModel)
    get(model.data, :infos, Dict())
end

function get_options(model::ASModel; options=Dict())
    opts = get(model.data, :options, Dict())
    return rmerge(opts, options)
end

function get_definitions(model::ASModel)
    # parse defs so values are Expr
    defs = get(model.data,:definitions,Dict())
    dynvars = cat(1, get_variables(model), keys(defs)...)
    _defs = OrderedDict{Symbol,Expr}([(Symbol(k), Dolo.sanitize( Dolo._to_expr(v),dynvars)) for (k, v) in defs])
    return _defs
end

function get_calibration(model::ASModel)
    calib = get(model.data, :calibration, Dict())
    # prep calib: parse to Expr, Symbol, or Number
    _calib  = OrderedDict{Symbol,Union{Expr,Symbol,Number}}()
    # add calibration for definitions
    for k in keys(calib)
        _calib[Symbol(k)] = Dolo._expr_or_number(calib[k])
    end
    # so far _calib is a symbolic calibration
    calibration = solve_triangular_system(_calib) ::OrderedDict{Symbol,Real }
    symbols = get_symbols( model )
    return ModelCalibration( calibration, symbols )
end

function get_domain(model::ASModel)
    domain = deepcopy(get(model.data, :domain, Dict()))
    calib = get_calibration(model)
    # TODO deal with empty dict and make robust construction
    states = get_symbols(model)[:states]
    min = [eval_with(calib, domain[(k)][1]) for k in states]
    max = [eval_with(calib, domain[(k)][2]) for k in states]
    return Domain(states, min, max)
end

function get_grid(model::ASModel; options=Dict())
    domain = get_domain(model)
    d = length(domain.states)
    options = get_options(model; options=Dict(:grid=>options))
    grid_dict = options[:grid]
    if :type in keys(grid_dict)
        grid_dict[:tag] = grid_dict[:type]
    end
    if grid_dict[:tag] == :Cartesian
        orders = get(grid_dict, :orders, [20 for i=1:d])
        grid = Dolo.CartesianGrid(domain.min, domain.max, orders)
    elseif grid_dict[:tag] == :Smolyak
        mu = grid_dict[:mu]
        grid = Dolo.SmolyakGrid(domain.min, domain.max, mu)
    else
        error("Unknown grid type.")
    end
    return grid
end

function get_exogenous(model::ASModel)
    exo_dict = get(model.data,:exogenous,Dict{Symbol,Any}())
    if length(exo_dict)==0
        exo_dict = get(model.data[:options], :exogenous, Dict{Symbol,Any}())
    end
    exo = Dolo._symbol_dict(exo_dict)
    calib = get_calibration(model)
    exogenous = Dolo._build_exogenous_entry(exo, calib)
    return exogenous
end

function set_calibration(model::ASModel, key::Symbol, value::Union{Real,Expr, Symbol})
    model.data[:calibration][key] = value
end


type Model{ID}<:AModel{ID}

    data
    name::String           # weakly immutable
    filename::String
    symbols::OrderedDict{Symbol,Array{Symbol,1}}        # weakly immutable
    equations      # weakly immutable
    definitions    # weakly immutable
    # functions      # weakly immutable
    factories::Dict
    calibration::ModelCalibration
    exogenous
    domain
    grid
    options

    function Model(data; print_code=false, filename="<string>")

        model = new(data)
        model.name = get_name(model)
        model.filename = filename
        model.symbols = get_symbols(model)
        model.equations = get_equations(model)
        model.definitions = get_definitions(model)
        model.factories = Dict()
        model.calibration = get_calibration(model)
        model.exogenous = get_exogenous(model)
        model.domain = get_domain(model)
        model.options = get_options(model)
        model.grid = get_grid(model)

        # now let's compile the functions:
        for eqtype in keys(model.equations)
            factory = Dolang.FunctionFactory(model,eqtype)
            model.factories[eqtype] = factory
            code = make_method(factory)
            print_code && println(code)
            eval(Dolo, code)
        end

        return model
    end

end

_numeric_mod_type{ID}(::Model{ID}) = Model{ID}


function Base.show(io::IO, model::Model)
    @printf("Model \n")
end

#### import functions


function Model(url::AbstractString; print_code=false)
    # it looks like it would be cool to use the super constructor ;-)
    if match(r"(http|https):.*", url) != nothing
        res = get(url)
        buf = Dolo.IOBuffer(res.data)
        data = Dolo._symbol_dict(Dolo.load(buf, yaml_types))
    else
        data = Dolo._symbol_dict(Dolo.load_file(url, yaml_types))
    end
    id = gensym()
    fname = basename(url)
    return Model{id}(data; print_code=print_code, filename=fname)
end


"""
Imports the model from a yaml file specified by the `url` input
parameter, and returns the corresponding `Dolo` model object.
"""
function yaml_import(url; print_code::Bool=false)
    Model(url; print_code=print_code)
end
