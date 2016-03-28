# -------------- #
# Symbolic model #
# -------------- #

immutable SymbolicModel <: ASM
    symbols::OrderedDict{Symbol,Vector{Symbol}}
    equations::OrderedDict{Symbol,Vector{Expr}}
    calibration::OrderedDict{Symbol,Union{Expr,Symbol,Number}}
    options::Dict{Symbol,Any}
    distribution::Dict{Symbol,Any}
    model_type::Symbol
    name::UTF8String
    filename::UTF8String

    function SymbolicModel(recipe::Associative, symbols::Associative,
                           eqs::Associative, calib::Associative,
                           options::Associative, dist::Associative,
                           name="modeldoesnotwork", filename="none")
        # prep symbols
        model_type = symbol(recipe[:model_spec])
        _symbols = OrderedDict{Symbol,Vector{Symbol}}()
        for k in recipe[:symbols]
            _symbols[symbol(k)] = Symbol[symbol(v) for v in get(symbols, k, [])]
        end

        # prep equations: parse to Expr
        _eqs = OrderedDict{Symbol,Vector{Expr}}()
        for k in keys(recipe[:specs])
            k == :arbitrage && continue  # we handle these separately

            these_eq = get(eqs, string(k), [])

            # verify that we have at least 1 equation if section is required
            if !get(recipe[:specs][k], :optional, false)
                length(these_eq) == 0 && error("equation section $k required")
            end

            # finally pass in the expressions
            _eqs[k] = Expr[_to_expr(eq) for eq in these_eq]
        end

        # handle the arbitrage, controls_lb, and controls_ub separately
        if haskey(recipe[:specs], :arbitrage)
            c_lb, c_ub, arb = _handle_arbitrage(eqs["arbitrage"],
                                                _symbols[:controls])
            _eqs[:arbitrage] = arb
            _eqs[:controls_lb] = c_lb
            _eqs[:controls_ub] = c_ub
        end

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

        new(_symbols, _eqs, _calib, dist, options, model_type, name, filename)
    end
end

function Base.show(io::IO, sm::SymbolicModel)
    println(io, """SymbolicModel
    - name: $(sm.name)
    """)
end

function SymbolicModel(from_yaml::Dict, model_type::Symbol, filename="none")
    # verify that we have all the required fields
    for k in ("symbols", "equations", "calibration")
        if !haskey(from_yaml, k)
            error("Yaml file must define section $k for DTCSCC model")
        end
    end

    d = deepcopy(from_yaml)
    recipe = RECIPES[model_type]
    out = SymbolicModel(recipe, pop!(d, "symbols"),
                        pop!(d, "equations"),
                        pop!(d, "calibration"),
                        _symbol_dict(pop!(d, "options", Dict())),
                        _symbol_dict(pop!(d, "distribution", Dict())),
                        pop!(d, "name", "modeldoesnotwork"),
                        filename)

    if !isempty(d)
        m = string("Fields $(join(keys(d), ", ", ", and ")) from yaml file ",
                   " were not used when constructing SymbolicModel")
        warn(m)
    end
    out
end

# add this method to be consistent with `model_type(::ANM)`
model_spec(sm::SymbolicModel) = sm.model_type

# ----------- #
# Calibration #
# ----------- #

# wrapper around ordered dict that will allow us to implement our own
# getindex/setindex! methods
immutable ModelCalibration
    flat::OrderedDict{Symbol,Float64}
    grouped::Dict{Symbol,Vector{Float64}}
end

function ModelCalibration(sm::SymbolicModel)
    flat = solve_triangular_system(sm)
    grouped = Dict{Symbol,Vector{Float64}}()
    for (k, nms) in sm.symbols
        grouped[k] = [flat[nm] for nm in nms]
    end

    ModelCalibration(flat, grouped)
end

# TODO: Decide if we should keep these semantics. Right now I've implemented
#       it so that if we do mc[x::Symbol] we try to extract the calibrated
#       value for a single parameter as a Float64. Then if we call
#       mc[s::AbstractString] we try to extract a whole symbol group as a
#       Vector{Float64}
Base.getindex(mc::ModelCalibration, n::Symbol) = mc.flat[n]
Base.getindex(mc::ModelCalibration, n::AbstractString) = mc.grouped[symbol(n)]

# now define methods that let us extract multiple params or groups at a time
Base.getindex(mc::ModelCalibration, n1::Symbol, nms::Symbol...) =
    [mc[n] for n in vcat(n1, nms...)]
Base.getindex(mc::ModelCalibration, n1::AbstractString, nms::AbstractString...) =
    Vector{Float64}[mc[n] for n in vcat(n1, nms...)]

# -------------------- #
# Model specific types #
# -------------------- #
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

# NOTE: a type parameter is needed so the `functions` field is not abstract
#       when instances of this type are actually created.
immutable DTCSCCModel{_T<:DTCSCCfunctions} <: ANM
    symbolic::SymbolicModel
    functions::_T
    calibration::ModelCalibration
end

immutable DTMSCCfunctions{T1,T2,T3,T4,T5,T6,T7,T8,T9}
    arbitrage::T1
    transition::T2
    auxiliary::T3
    value::T4
    expectation::T5
    direct_response::T6
    controls_lb::T7
    controls_ub::T8
    arbitrage_2::T9
end

immutable DTMSCCModel{_T<:DTMSCCfunctions} <: ANM
    symbolic::SymbolicModel
    functions::_T
    calibration::ModelCalibration
end

for (TF, TM, ms) in [(:DTCSCCfunctions, :DTCSCCModel, :(:dtcscc)),
                     (:DTMSCCfunctions, :DTMSCCModel, :(:dtmscc))]
    @eval begin
        model_spec(::$(TM)) = $ms
        model_spec(::Type{$(TM)}) = $ms
        model_spec(::$(TF)) = $ms
        model_spec(::Type{$(TF)}) = $ms

        # function type constructor
        function $(TF)(sm::SymbolicModel)
            if model_spec(sm) != model_spec($TF)
                msg = string("Symbolic model is of type $(model_spec(sm)) ",
                             "cannot create functions of type $($TF)")
                error(msg)
            end
            $(TF)([let
                       eval(compile_equation(sm, fld))
                   end
                   for fld in fieldnames($(TF))]...)
        end

        # model type constructor
        function Base.convert(::Type{$TM}, sm::SymbolicModel)
            if model_spec(sm) != model_spec($TM)
                msg = string("Symbolic model is of type $(model_spec(sm)) ",
                             "cannot create model of type $($TM)")
                error(msg)
            end
            $(TM)(sm, $(TF)(sm), ModelCalibration(sm))
        end

        Base.convert(::Type{SymbolicModel}, m::$(TM)) = m.symbolic
    end
end
