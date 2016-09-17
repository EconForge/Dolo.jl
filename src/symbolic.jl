# -------------- #
# Symbolic Model #
# -------------- #

immutable SymbolicModel{ID} <: ASM{ID}
    symbols::OrderedDict{Symbol,Vector{Symbol}}
    equations::OrderedDict{Symbol,Vector{Expr}}
    calibration::OrderedDict{Symbol,Union{Expr,Symbol,Number}}
    exogenous::Dict{Symbol,Any}
    options::Dict{Symbol,Any}
    definitions::OrderedDict{Symbol,Expr}
    name::String
    filename::String

    function SymbolicModel(recipe::Associative, symbols::Associative,
                           eqs::Associative, calib::Associative,
                           exog::Associative,
                           options::Associative, defs::Associative,
                           name="modeldoesnotwork", filename="none")
        # prep symbols
        model_type = Symbol(recipe[:model_spec])
        _symbols = OrderedDict{Symbol,Vector{Symbol}}()
        for _ in recipe[:symbols]
            k = Symbol(_)
            _symbols[k] = Symbol[Symbol(v) for v in get(symbols, k, [])]
        end

        # prep equations: parse to Expr
        _eqs = OrderedDict{Symbol,Vector{Expr}}()
        for k in keys(recipe[:specs])

            # we handle these separately
            (k in [:arbitrage,]) && continue

            these_eq = get(eqs, k, [])

            # verify that we have at least 1 equation if section is required
            if !get(recipe[:specs][k], :optional, false)
                length(these_eq) == 0 && error("equation section $k required")
            end

            # finally pass in the expressions
            _eqs[k] = Expr[_to_expr(eq) for eq in these_eq]
        end

        # handle the arbitrage, arbitrage_exp, controls_lb, and controls_ub
        if haskey(recipe[:specs], :arbitrage)
            c_lb, c_ub, arb = _handle_arbitrage(eqs[:arbitrage],
                                                _symbols[:controls])
            _eqs[:arbitrage] = arb
            _eqs[:controls_lb] = c_lb
            _eqs[:controls_ub] = c_ub
        end

        # parse defs so values are Expr
        _defs = OrderedDict{Symbol,Expr}([(k, _to_expr(v)) for (k, v) in defs])

        # prep calib: parse to Expr, Symbol, or Number
        _calib  = OrderedDict{Symbol,Union{Expr,Symbol,Number}}()
        for k in keys(_symbols)
            for nm in _symbols[k]
                if k == :shocks
                    _calib[nm] = 0.0
                else
                    _calib[nm] = _expr_or_number(calib[nm])
                end
            end
        end

        # add calibration for definitions
        for k in keys(_defs)
            _calib[k] = _expr_or_number(calib[k])
        end

        new(_symbols, _eqs, _calib, exog, options, _defs, name, filename)
    end
end

function Base.show(io::IO, sm::SymbolicModel)
    println(io, """SymbolicModel
    - name: $(sm.name)
    """)
end

function SymbolicModel(data::Dict, filename="none")
    # verify that we have all the required fields
    for k in (:symbols, :equations, :calibration)
        if !haskey(data, k)
            error("Yaml file must define section $k for dtcc model")
        end
    end

    d = _symbol_dict(deepcopy(data))
    mt = pop!(d, :model_type, :dtcc)
    Symbol(mt) == :dtcc || error("Only support dtcc models now")

    recipe = RECIPES[:dtcc]
    nm = pop!(d, :name, "modeldoesnotwork")
    id = gensym(nm)
    options = pop!(d, :options, Dict{Symbol,Any}())
    defs = pop!(d, :definitions, Dict{Symbol,Any}())
    # exog = pop!(d, :exogenous, Dict{Symbol,Any}())
    exog = get(options, :exogenous, Dict{Symbol,Any}())
    if exog == Dict{Symbol,Any}()
      error("Yaml file must define section 'exogenous' for dtcc model")
    end
    out = SymbolicModel{id}(recipe, pop!(d, :symbols),
                            pop!(d, :equations),
                            pop!(d, :calibration),
                            exog,
                            options,
                            defs,
                            nm,
                            filename)

    if !isempty(d)
        m = string("Fields $(join(keys(d), ", ", ", and ")) from yaml file ",
                   " were not used when constructing SymbolicModel")
        warn(m)
    end
    out
end

# assume dolo style model is default, special case others
function _get_args(sm::SymbolicModel, spec)
    # get args
    args = OrderedDict{Symbol,Vector{Tuple{Symbol,Int}}}()
    for (grp, shift, nm) in spec[:eqs]
        grp == "parameters" && continue

        syms = sm.symbols[Symbol(grp)]
        args[Symbol(nm)] = Tuple{Symbol,Int}[(v, shift) for v in syms]
    end
    args
end

function Dolang.FunctionFactory(sm::SymbolicModel, func_nm::Symbol)
    spec = RECIPES[:dtcc][:specs][func_nm]
    eqs = sm.equations[func_nm]

    # get targets
    target = get(spec, :target, [nothing])[1]
    has_targets = !(target === nothing)
    targets = has_targets ? sm.symbols[Symbol(target)] : Symbol[]

    # get other stuff
    args = Dolo._get_args(sm, spec)
    params = sm.symbols[:parameters]
    dispatch = _numeric_mod_type(sm)

    FunctionFactory(dispatch, eqs, args, params, targets=targets,
                    defs=sm.definitions, funname=func_nm)
end

# ----- #
# Tools #
# ----- #

# Used for extracting complementarity conditions from arbitrage equations
function _handle_arbitrage(arb, controls)
    controls_lb = Expr[]
    controls_ub = Expr[]
    arbitrage = Expr[]
    for (i, v) in enumerate(arb)
        parts = split(v, "|")

        if length(parts) == 1
            push!(arbitrage, _to_expr(parts[1]))
            push!(controls_lb, _to_expr(-Inf))
            push!(controls_ub, _to_expr(Inf))

        # We have complementarity conditions to deal with
        elseif length(parts) == 2
            # convert arbitrage equation to an expression
            push!(arbitrage, _to_expr(parts[1]))

            # Now deal with complementarities
            c_parts = split(parts[2], "<=")
            n_c = length(c_parts)

            if n_c == 3
                push!(controls_lb, _to_expr(inf_to_Inf(parse(c_parts[1]))))
                push!(controls_ub, _to_expr(inf_to_Inf(parse(c_parts[3]))))

                # verify that the control is what we want
                mid = parse(c_parts[2])
                if mid != controls[i]
                    msg = string("Error in complementarity condition. ",
                                 "Expected $(controls[i]) found $mid")
                    error(msg)
                end
            elseif n_c == 2
                # only have one_sided condition. Need to work a bit harder
                ind = 0
                for (j, ex) in enumerate(c_parts)
                    if parse(ex) == controls[i]
                        ind = j
                        break
                    end
                end

                # we didn't find the control, throw an error
                if ind == 0
                    msg = string("Error in complementarity condition. ",
                                 "Expected $(controls[i]), but did not find.")
                    error(msg)
                end

                # otherwise continue on
                if ind == 1
                    push!(controls_lb, _to_expr(-Inf))
                    push!(controls_ub, _to_expr(inf_to_Inf(parse(c_parts[2]))))
                else
                    push!(controls_lb, _to_expr(inf_to_Inf(parse(c_parts[1]))))
                    push!(controls_ub, _to_expr(Inf))
                end
            else
                msg = string("Malformed complementarity condition. ",
                             "Expected 1 or 2 `<=`, found $(n_c-1)")
                error(msg)
            end

        else

            msg = string("Malformed arbitrage equation. ",
                         "Only have one occurance of `|` allowed.")
            error(msg)
        end

    end
    controls_lb, controls_ub, arbitrage
end
