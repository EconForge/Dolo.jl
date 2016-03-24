# --------------- #
# Utility methods #
# --------------- #

_to_expr(x::Expr) = x
_to_expr(x::Union{Symbol,Number}) = Expr(:block, x)
_to_expr(x::AbstractString) = parse(x)

_expr_or_number(x::Union{AbstractString,Symbol,Expr}) = _to_expr(x)
_expr_or_number(x::Number) = x

inf_to_Inf(x::Number) = x
inf_to_Inf(x::Symbol) = @match x begin
    inf => Inf
    -inf => -Inf
    _ => x
end

function solve_triangular_system(sm::ASM)
    solutions = Dict{Symbol,Number}()
    finished = false
    dict = sm.calibration
    N = length(dict)
    n = 0

    while !(finished || n>N)
        done_smthg = false
        n += 1
        for k in keys(dict)
            if !haskey(solutions, k)
                expr = dict[k]
                try
                    sol = eval(:(let $([:($x=$y) for (x, y) in solutions]...); $expr end))
                    solutions[k] = sol
                    done_smthg = true
                end
            end
        end
        if done_smthg == false
            finished = true
        end
    end

    if length(solutions) < length(dict)
        error("Not a triangular system")
    end

    # reorder solutions to match sm.calibration
    OrderedDict{Symbol,Number}([(k, solutions[k]) for k in keys(sm.calibration)])
end

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
                             "Expected 2 `<=`, found $(n_c-1)")
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

# ------------------------------ #
# Symboloc model and constructor #
# ------------------------------ #

immutable DTCSCCSymbolicModel <: ASM
    symbols::OrderedDict{Symbol,Vector{Symbol}}
    equations::OrderedDict{Symbol,Vector{Expr}}
    calibration::OrderedDict{Symbol,Union{Expr,Symbol,Number}}
    options::Dict{Symbol,Any}
    distribution::Dict{Symbol,Any}
    name::UTF8String
    filename::UTF8String

    function DTCSCCSymbolicModel(symbols::Associative, eqs::Associative,
                                 calib::Associative, options::Associative,
                                 dist::Associative,
                                 name="modeldoesnotwork", filename="none")

        # prep symbols
        _symbols = OrderedDict{Symbol,Vector{Symbol}}()
        for k in map(symbol, RECIPES[:dtcscc][:symbols])
            _symbols[k] = Symbol[symbol(v) for v in get(symbols, string(k), [])]
        end

        # prep equations: parse to Expr
        _eqs = OrderedDict{Symbol,Vector{Expr}}()
        for k in map(symbol, keys(RECIPES[:dtcscc][:specs]))
            k == :arbitrage && continue  # we handle these separately

            these_eq = get(eqs, string(k), [])

            # verify that we have at least 1 equation if section is required
            if !get(RECIPES[:dtcscc][:specs][k], :optional, false)
                length(these_eq) == 0 && error("equation section $k required")
            end

            # finally pass in the expressions
            _eqs[k] = Expr[_to_expr(eq) for eq in these_eq]
        end

        # handle the arbitrage, controls_lb, and controls_ub separately
        c_lb, c_ub, arb = _handle_arbitrage(eqs["arbitrage"],
                                            _symbols[:controls])
        _eqs[:arbitrage] = arb
        _eqs[:controls_lb] = c_lb
        _eqs[:controls_ub] = c_ub

        # prep calib: parse to Expr, Symbol, or Number
        _calib  = OrderedDict{Symbol,Union{Expr,Symbol,Number}}()
        for k in keys(_symbols)
            for nm in _symbols[k]
                if k == :shocks
                    _calib[nm] = 0.0
                else
                    _calib[nm] = _expr_or_number(calib[string(nm)])
                end
            end
        end

        _distribution

        new(_symbols, _eqs, _calib, dist, options, name, filename)
    end
end

function DTCSCCSymbolicModel(from_yaml::Dict, filename="none")
    # verify that we have all the required fields
    for k in ("symbols", "equations", "calibration")
        if !haskey(from_yaml, k)
            error("Yaml file must define section $k for DTCSCC model")
        end
    end

    d = copy(from_yaml)
    out = DTCSCCSymbolicModel(pop!(d, "symbols"),
                              pop!(d, "equations"),
                              pop!(d, "calibration"),
                              _symbol_dict(pop!(d, "options", Dict())),
                              _symbol_dict(pop!(d, "distribution", Dict())),
                              pop!(d, "name", "modeldoesnotwork"),
                              filename)

    !isempty(d) && warn("Fields $(collect(keys(d))) from yaml file were not used")
    out
end

# ----------------------------- #
# Numeric Model and constructor #
# ----------------------------- #
# TODO: decide if we really want all 8 type params. It doesn't hurt, it just
#       looks funny
immutable DTCSCCfunctions{T1,T2,T3,T4,T5,T6,T7,T8}
    arbitrage::T1
    transition::T2
    auxiliary::T3
    value::T4
    expectation::T5
    direct_response::T6
    controls_lb::T7
    controls_ub::T8
end

DTCSCCfunctions(sm::DTCSCCSymbolicModel) =
    DTCSCCfunctions([eval(compile_equation(sm, fld))
                     for fld in fieldnames(DTCSCCfunctions)]...)

# wrapper around ordered dict that will allow us to implement our own
# getindex/setindex! methods
immutable ModelCalibration
    data::OrderedDict{Symbol,Number}
end

ModelCalibration(sm::DTCSCCSymbolicModel) =
    ModelCalibration(solve_triangular_system(sm))

# NOTE: a type parameter is needed so the `functions` field is not abstract
#       when instances of this type are actually created.
immutable DTCSCCModel{_T<:DTCSCCfunctions} <: ANM
    symbolic::DTCSCCSymbolicModel
    functions::_T
    calibration::ModelCalibration
end

Base.convert(::Type{DTCSCCModel}, sm::DTCSCCSymbolicModel) =
    DTCSCCModel(sm, DTCSCCfunctions(sm), ModelCalibration(sm))

Base.convert(::Type{DTCSCCSymbolicModel}, m::DTCSCCModel) = m.symbolic

for T in (:DTCSCCModel, :DTCSCCSymbolicModel)
    @eval begin
        model_spec(::$(T)) = :dtcscc
        model_spec(::Type{$(T)}) = :dtcscc
    end
end
