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

        new(_symbols, _eqs, _calib, options, dist, model_type, name, filename)
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

# ----------- #
# Calibration #
# ----------- #

immutable FlatCalibration <: Associative{Symbol,Float64}
    d::OrderedDict{Symbol,Float64}
end

FlatCalibration(pairs::Pair{Symbol,Float64}...) =
    FlatCalibration(OrderedDict(pairs))

FlatCalibration() = FlatCalibration(OrderedDict{Symbol,Float64}())

immutable GroupedCalibration <: Associative{Symbol,Vector{Float64}}
    d::Dict{Symbol,Vector{Float64}}
end

GroupedCalibration() = GroupedCalibration(Dict{Symbol,Vector{Float64}}())

GroupedCalibration(pairs::Pair{Symbol,Vector{Float64}}...) =
    GroupedCalibration(Dict(pairs))

# use let block so `Calib` isn't a member of the Dolo module
let
    typealias Calib Union{FlatCalibration,GroupedCalibration}

    Base.getindex(c::Calib, n::Symbol) = c.d[n]
    Base.getindex(c::Calib, nms::Symbol...) = [c[n] for n in nms]

    # define the rest of the Associative interface by forwarding to c.d

    # 1-arg functions
    for f in [:length, :keys, :values, :keytype, :valtype, :eltype, :isempty,
              :start]
        @eval Base.$(f)(c::$(Calib)) = $(f)(c.d)
    end

    # 2-arg functions
    for f in [:haskey, :delete!, :pop!, :sizehint!, :next, :done]
        @eval Base.$(f)(c::$(Calib), arg) = $(f)(c.d, arg)
    end

    # 3-arg functions
    for f in [:get, :get!, :getkey, :pop!]
        @eval Base.$(f)(c::$(Calib), arg1, arg2) = $(f)(c.d, arg1, arg2)
    end

    # methods that didn't match signatures above
    Base.get!(f::Function, c::Calib, key) = get!(f, c.d, key)
    Base.get(f::Function, c::Calib, key) = get(f, c.d, key)
    Base.merge!{T<:Calib}(c::T, others::T...) =
        T(merge!(c.d, [o.d for o in others]...))
    Base.show(io::IO, c::Calib) = show(io, c.d)
    Base.showdict(io::IO, c::Calib) = Base.showdict(io, c.d)
    Base.copy{T<:Calib}(c::T) = T(copy(c.d))
    Base.deepcopy{T<:Calib}(c::T) = T(deepcopy(c.d))
end

# setting a single value
Base.setindex!(fc::FlatCalibration, v::Real, k::Symbol) =
    fc.d[k] = convert(Float64, v)

# setting multiple values
function Base.setindex!(fc::FlatCalibration, vs::Union{Tuple,AbstractVector},
                        ks::Symbol...)
    if length(vs) != length(ks)
        throw(DimensionMismatch("length of keys and values must be the same"))
    end

    for (v, k) in zip(vs, ks)
        fc[k] = v
    end
    fc
end

# setting a single group
function Base.setindex!(gc::GroupedCalibration, v::AbstractVector, k::Symbol)
    # use get with default of v in case we are adding a new key
    if length(v) != length(get(gc, k, v))
        msg = string("length $(length(v)) of new values for key $k does",
                     "not match length of existing value ($(length(gc[k])))")
        throw(DimensionMismatch(msg))
    end
    gc.d[k] = convert(Vector{Float64}, v)
end

# setting multiple groups
function Base.setindex!(gc::GroupedCalibration, vs::Tuple, ks::Symbol...)
    if length(vs) != length(ks)
        throw(DimensionMismatch("length of keys and values must be the same"))
    end

    for (v, k) in zip(vs, ks)
        gc[k] = v
    end
    gc
end

# Contains a FlatCalibration and GroupedCalibration, plus information on symbols
immutable ModelCalibration
    flat::FlatCalibration
    grouped::GroupedCalibration
    symbol_table::Dict{Symbol,Tuple{Symbol,Int}}
    symbol_groups::OrderedDict{Symbol,Vector{Symbol}}
end

function ModelCalibration(sm::SymbolicModel)
    flat = FlatCalibration(solve_triangular_system(sm))
    grouped = GroupedCalibration()
    for (k, nms) in sm.symbols
        grouped[k] = Float64[flat[nm] for nm in nms]
    end

    symbol_table = Dict{Symbol,Tuple{Symbol,Int}}()
    for (grp, vals) in sm.symbols
        for (i, v) in enumerate(vals)
            symbol_table[v] = (grp, i)
        end
    end

    # make sure we documented where in grouped every symbol is
    @assert sort(collect(keys(symbol_table))) == sort(collect(keys(flat)))

    ModelCalibration(flat, grouped, symbol_table, deepcopy(sm.symbols))
end

for f in (:copy, :deepcopy)
    @eval Base.$(f)(mc::ModelCalibration) =
        ModelCalibration($(f)(mc.flat), $(f)(mc.grouped),
                         $(f)(mc.symbol_table), $(f)(mc.symbol_groups))
end

# getindex and setindex! should be forwarded on to the grouped field
Base.getindex(mc::ModelCalibration, n::Symbol) = mc.grouped[n]
Base.getindex(mc::ModelCalibration, nms::Symbol...) = [mc[n] for n in nms]

function _setindex_group!(mc::ModelCalibration, v::AbstractVector, grp::Symbol)
    ks = mc.symbol_groups[grp]
    if length(v) != length(ks)
        msg = string("Calibration has $(length(ks)) symbols in $grp, ",
        "but passed $(length(v)) values")
        throw(DimensionMismatch(msg))
    end

    # update grouped
    mc.grouped[grp] = v

    # update flat
    for (vi, ki) in zip(v, ks)
        mc.flat[ki] = vi
    end
    mc
end

function _setindex_flat!(mc::ModelCalibration, v::Real, k::Symbol)
    # update flat
    mc.flat[k] = v

    # update grouped
    grp, ix = mc.symbol_table[k]
    mc.grouped[grp][ix] = v
    mc
end

function Base.setindex!(mc::ModelCalibration, v, k::Symbol)
    if haskey(mc.symbol_table, k)       # setting a variable
        _setindex_flat!(mc, v, k)
    elseif haskey(mc.symbol_groups, k)  # setting a group
        _setindex_group!(mc, v, k)
    else
        throw(KeyError("ModelCalibration has no variable or group named $k"))
    end
end

function Base.setindex!(mc::ModelCalibration, vs::Tuple, ks::Symbol...)
    if length(vs) != length(ks)
        throw(DimensionMismatch("length of keys and values must be the same"))
    end

    for (k, v) in zip(ks, vs)
        mc[k] = v
    end
    mc
end

# tries to replace a symbol if the key is in the calibration, otherwise just
# keeps the symbol in place
_replace_me(mc::ModelCalibration, s::Symbol) = get(mc.flat, s, s)
_replace_me(mc, o) = o

# eval with will work on
function eval_with(mc::ModelCalibration, ex::Expr)
    # put in let block to allow us to define intermediates in expr and not
    # have them become globals in `current_module()` at callsite
    new_ex = MacroTools.prewalk(s->_replace_me(mc, s), ex)
    eval(Dolo, :(
    let
        $new_ex
    end))
end

eval_with(mc::ModelCalibration, s::AbstractString) = eval_with(mc, _to_expr(s))
eval_with(mc::ModelCalibration, s::Symbol) = _replace_me(mc, s)
eval_with(mc::ModelCalibration, d::Associative) =
    Dict{Symbol,Any}([(symbol(k), eval_with(mc, v)) for (k, v) in d])
eval_with(mc::ModelCalibration, x::Number) = x
eval_with(mc::ModelCalibration, x::AbstractArray) = map(y->eval_with(mc, y), x)

# ------------------- #
# Numeric model types #
# ------------------- #
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
    options::Dict{Symbol,Any}
    distribution::Dict{Symbol,Any}
    model_type::Symbol
    name::UTF8String
    filename::UTF8String
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
    options::Dict{Symbol,Any}
    distribution::Dict{Symbol,Any}
    model_type::Symbol
    name::UTF8String
    filename::UTF8String
end

for (TF, TM, ms) in [(:DTCSCCfunctions, :DTCSCCModel, :(:dtcscc)),
                     (:DTMSCCfunctions, :DTMSCCModel, :(:dtmscc))]
    @eval begin
        model_type(::$(TM)) = $ms
        model_type(::Type{$(TM)}) = $ms
        model_type(::$(TF)) = $ms
        model_type(::Type{$(TF)}) = $ms

        # function type constructor
        function $(TF)(sm::SymbolicModel; print_code::Bool=false)
            if model_type(sm) != model_type($TF)
                msg = string("Symbolic model is of type $(model_type(sm)) ",
                             "cannot create functions of type $($TF)")
                error(msg)
            end
            $(TF)([let
                       eval(Dolo, compile_equation(sm, fld; print_code=print_code))
                   end
                   for fld in fieldnames($(TF))]...)
        end
        Base.convert(::Type{SymbolicModel}, m::$(TM)) = m.symbolic

        # model type constructor
        function ($TM)(sm::SymbolicModel; print_code::Bool=false)
            if model_type(sm) != model_type($TM)
                msg = string("Symbolic model is of type $(model_type(sm)) ",
                             "cannot create model of type $($TM)")
                error(msg)
            end
            calib = ModelCalibration(sm)
            options = eval_with(calib, deepcopy(sm.options))
            dist = eval_with(calib, deepcopy(sm.distribution))
            # hack to parse normal transition matrix into a matrix instead of
            # a vector of vectors
            if haskey(dist, :Normal)
                n = length(calib[:shocks])
                dist[:Normal] = reshape(vcat(dist[:Normal]...), n, n)
            end
            funcs = $(TF)(sm; print_code=print_code)
            $(TM)(sm, funcs, calib, options, dist, sm.model_type,
                  sm.name, sm.filename)
        end
    end
end

# ------------- #
# Other methods #
# ------------- #

model_type(m::AbstractModel) = m.model_type
filename(m::AbstractModel) = m.filename
name(m::AbstractModel) = m.name
