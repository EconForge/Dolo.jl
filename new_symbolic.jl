import Dolo
using YAML
url = "examples/models/rbc_dtcc_mc.yaml"


# used for constructing appropraite dict from YAML object.
function construct_type_map(t::Symbol, constructor::YAML.Constructor,
                            node::YAML.Node)
    return node
end
pairs = [("!Cartesian", :Cartesian),
         ("!Normal", :Normal),
         ("!MarkovChain", :MarkovChain),
         ("!AR1", :AR1)]
yaml_types = Dict{AbstractString,Function}([(t, (c, n) -> construct_type_map(s, c, n))
                       for (t, s) in pairs])
data = Dolo.load_file(url, yaml_types)

data
using DataStructures

module NewDolo

    import Dolo
    using DataStructures

    type Model
        data::Dict{Any,Any}
    end

    function get_symbols(model::Model)
        # model_type = Symbol(recipe[:model_spec])
        # _symbols = OrderedDict{Symbol,Vector{Symbol}}()
        # for _ in recipe[:symbols]
        #     k = Symbol(_)
        #     _symbols[k] = Symbol[Symbol(v) for v in get(symbols, k, [])]
        # end
        syms = model.data["symbols"]
        symbols = Dict()
        for k in keys(syms)
            symbols[Symbol(k)] = [Symbol(e) for e in syms[k]]
        end
        return symbols
    end

    function get_name(model::Model)
        get(model.data, "name", "modeldoesnotwork")
    end

    function get_equations(model::Model)

        eqs = model.data["equations"]
        println(eqs)
        recipe = Dolo.RECIPES[:dtcc]

        # prep equations: parse to Expr
        _eqs = OrderedDict{Symbol,Vector{Expr}}()
        for k in keys(recipe[:specs])
            # we handle these separately
            (k in [:arbitrage,]) && continue
            these_eq = get(eqs, string(k), [])
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
        model.data["definitions"]
    end

    function get_infos(model::Model)
        get(model.data, "infos", Dict())
    end

    function get_options(model::Model)
        get(model.data, "options", Dict())
    end

    function get_definitions(model::Model)
        # parse defs so values are Expr
        defs = get(model.data,"definitions",Dict())
        println(defs)
        _defs = OrderedDict{Symbol,Expr}([(Symbol(k), Dolo._to_expr(v)) for (k, v) in defs])
        return _defs
    end

    function get_calibration(model::Model)
        calib = get(model.data, "calibration", Dict())
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


end


import NewDolo

model = NewDolo.Model(data)
NewDolo.get_name(model)
NewDolo.get_symbols(model)
NewDolo.get_definitions(model)
NewDolo.get_equations(model)
NewDolo.get_infos(model)
NewDolo.get_options(model)
NewDolo.get_definitions(model)
cal = NewDolo.get_calibration(model)
cal.flat
