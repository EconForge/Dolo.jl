import YAML
import DataStructures:OrderedDict


function construct_type_map(t::Symbol, constructor::YAML.Constructor,
                            node::YAML.Node)
    mapping = _symbol_dict(YAML.construct_mapping(constructor, node))
    mapping[:tag] = t
    mapping
end


pairs = [("!Cartesian", :Cartesian),
         ("!Smolyak", :Smolyak),
         ("!Normal", :Normal),
         ("!MarkovChain", :MarkovChain),
         ("!AR1", :AR1)]
yaml_types = Dict{AbstractString,Function}([(t, (c, n) -> construct_type_map(s, c, n))
                       for (t, s) in pairs])

abstract type AModel end
abstract type ANModel <: AModel end

type Model <:AModel
    data::Dict{Any,Any}
end

type Domain
    states::Vector{Symbol}
    min::Vector{Float64}
    max::Vector{Float64}
end

function Model(url::AbstractString)
    if match(r"(http|https):.*", url) != nothing
        res = get(url)
        buf = IOBuffer(res.data)
        data = _symbol_dict(load(buf, yaml_types))
    else
        data = _symbol_dict(load_file(url, yaml_types))
    end
    return Model(data)
end

function get_symbols(model::Model)
    syms = model.data[:symbols]
    symbols = OrderedDict{Symbol,Vector{Symbol}}()
    for k in keys(syms)
        symbols[Symbol(k)] = [Symbol(e) for e in syms[k]]
    end
    return symbols
end

function get_name(model::Model)
    get(model.data, :name, "modeldoesnotwork")
end

function get_equations(model::Model)

    eqs = model.data[:equations]
    recipe = RECIPES[:dtcc]

    # prep equations: parse to Expr
    _eqs = OrderedDict{Symbol,Vector{Expr}}()
    for k in keys(recipe[:specs])
        # we handle these separately
        (k in [:arbitrage,]) && continue
        these_eq = get(eqs, (k), [])
        # verify that we have at least 1 equation if section is required
        if !get(recipe[:specs][k], :optional, false)
            length(these_eq) == 0 && error("equation section $k required")
        end
        # finally pass in the expressions
        _eqs[k] = Expr[_to_expr(eq) for eq in these_eq]
    end
    # handle the arbitrage, arbitrage_exp, controls_lb, and controls_ub
    if haskey(recipe[:specs], "arbitrage")
        c_lb, c_ub, arb = _handle_arbitrage(eqs["arbitrage"],
                                            _symbols[:controls])
        _eqs[:arbitrage] = arb
        _eqs[:controls_lb] = c_lb
        _eqs[:controls_ub] = c_ub
    end

    return _eqs
end

function definitions(model::Model)
    model.data[:definitions]
end

function get_infos(model::Model)
    get(model.data, :infos, Dict())
end

function get_options(model::Model)
    get(model.data, :options, Dict())
end

function get_definitions(model::Model)
    # parse defs so values are Expr
    defs = get(model.data,:definitions,Dict())
    _defs = OrderedDict{Symbol,Expr}([(Symbol(k), _to_expr(v)) for (k, v) in defs])
    return _defs
end

function get_calibration(model::Model)
    calib = get(model.data, :calibration, Dict())
    # prep calib: parse to Expr, Symbol, or Number
    _calib  = OrderedDict{Symbol,Union{Expr,Symbol,Number}}()
    # add calibration for definitions
    for k in keys(calib)
        _calib[Symbol(k)] = _expr_or_number(calib[k])
    end
    # so far _calib is a symbolic calibration
    calibration = solve_triangular_system(_calib)
    symbols = get_symbols( model )
    println(typeof(calibration))
    println(typeof(symbols))
    return ModelCalibration( calibration, symbols )
end

function get_domain(model::Model)
    domain = deepcopy(get(model.data, :domain, Dict()))
    calib = get_calibration(model)
    # TODO deal with empty dict and make robust construction
    states = get_symbols(model)[:states]
    min = [eval_with(calib, domain[(k)][1]) for k in states]
    max = [eval_with(calib, domain[(k)][2]) for k in states]
    return Domain(states, min, max)
end

function get_grid(model::Model)
    domain = get_domain(model)
    d = length(domain.states)
    grid_dict = model.data[:options][:grid]
    if grid_dict[:tag] == :Cartesian
        # TODO simplify...

        ogrid = model.data[:options][:grid]
        orders = get(ogrid, :orders, [20 for i=1:d])
        grid = CartesianGrid(domain.min, domain.max, orders)
    elseif grid_dict[:tag] == :Smolyak
        ogrid = model.data[:options][:grid]
        mu = ogrid[:mu]
        grid = SmolyakGrid(domain.min, domain.max, mu)
    else
        error("Unknown grid type.")
    end

    return grid
end

function get_exogenous(model::Model)
    exo_dict = get(model.data,:exogenous,Dict{Symbol,Any}())
    if length(exo_dict)==0
        exo_dict = get(model.data[:options], :exogenous, Dict{Symbol,Any}())
    end
    exo = _symbol_dict(exo_dict)
    calib = get_calibration(model)
    exogenous = _build_exogenous_entry(exo, calib)
    return exogenous
end
