module NewDolo

    import YAML
    import Dolo
    using DataStructures


    function construct_type_map(t::Symbol, constructor::YAML.Constructor,
                                node::YAML.Node)
        mapping = Dolo._symbol_dict(YAML.construct_mapping(constructor, node))
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

    type Model
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
            buf = Dolo.IOBuffer(res.data)
            data = Dolo._symbol_dict(Dolo.load(buf, yaml_types))
        else
            data = Dolo._symbol_dict(Dolo.load_file(url, yaml_types))
        end
        return Model(data)
    end

    function get_symbols(model::Model)
        syms = model.data[:symbols]
        symbols = Dict()
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
        recipe = Dolo.RECIPES[:dtcc]

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
            _eqs[k] = Expr[Dolo._to_expr(eq) for eq in these_eq]
        end
        # handle the arbitrage, arbitrage_exp, controls_lb, and controls_ub
        if haskey(recipe[:specs], "arbitrage")
            c_lb, c_ub, arb = Dolo._handle_arbitrage(eqs["arbitrage"],
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
        _defs = OrderedDict{Symbol,Expr}([(Symbol(k), Dolo._to_expr(v)) for (k, v) in defs])
        return _defs
    end

    function get_calibration(model::Model)
        calib = get(model.data, :calibration, Dict())
        # prep calib: parse to Expr, Symbol, or Number
        _calib  = OrderedDict{Symbol,Union{Expr,Symbol,Number}}()
        # add calibration for definitions
        for k in keys(calib)
            _calib[Symbol(k)] = Dolo._expr_or_number(calib[k])
        end
        # so far _calib is a symbolic calibration
        calibration = Dolo.solve_triangular_system(_calib)
        symbols = get_symbols( model )
        return Dolo.ModelCalibration( calibration, symbols )
    end

    function get_domain(model::Model)
        domain = deepcopy(get(model.data, :domain, Dict()))
        calib = get_calibration(model)
        # TODO deal with empty dict and make robust construction
        states = get_symbols(model)[:states]
        min = [Dolo.eval_with(calib, domain[(k)][1]) for k in states]
        max = [Dolo.eval_with(calib, domain[(k)][2]) for k in states]
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
            grid = Dolo.CartesianGrid(domain.min, domain.max, orders)
        elseif grid_dict[:tag] == :Smolyak
            ogrid = model.data[:options][:grid]
            mu = ogrid[:mu]
            grid = Dolo.SmolyakGrid(domain.min, domain.max, mu)
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
        exo = Dolo._symbol_dict(exo_dict)
        calib = get_calibration(model)
        exogenous = Dolo._build_exogenous_entry(exo, calib)
        return exogenous
    end


end

import Dolo


import NewDolo


url = "examples/models/rbc_dtcc_mc.yaml"
typeof(url)

model = NewDolo.Model(url)
data = model.data


@time NewDolo.get_name(model)
@time NewDolo.get_symbols(model)
@time NewDolo.get_definitions(model)
@time NewDolo.get_equations(model)
@time NewDolo.get_infos(model)
@time NewDolo.get_options(model)
@time NewDolo.get_definitions(model)

# apparently the time is dominated by the solution of the calibration
# the result should be cached

@time calib = NewDolo.get_calibration(model)

@time domain = NewDolo.get_domain(model)
@time exo = NewDolo.get_exogenous(model)
@time grid = NewDolo.get_grid(model)
