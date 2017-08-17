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
    const Calib = Union{FlatCalibration,GroupedCalibration}

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

function ModelCalibration(calib::OrderedDict{Symbol,Real}, symbols::OrderedDict{Symbol,Array{Symbol,1}})
    flat = FlatCalibration(calib)
    grouped = GroupedCalibration()
    for (k, nms) in symbols
        grouped[k] = Float64[flat[nm] for nm in nms]
    end
    # grouped[:definitions] = Float64[flat[nm] for nm in keys(sm.definitions)]
    symbol_table = Dict{Symbol,Tuple{Symbol,Int}}()
    for (grp, vals) in symbols
        for (i, v) in enumerate(vals)
            symbol_table[v] = (grp, i)
        end
    end

    # for (i, v) in enumerate(keys(sm.definitions))
    #     symbol_table[v] = (:definitions, i)
    # end

    # # make sure we documented where in grouped every symbol is
    # @assert sort(collect(keys(symbol_table))) == sort(collect(keys(flat)))

    ModelCalibration(flat, grouped, symbol_table, deepcopy(symbols))
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
    # have them become globals in `Dolo`
    new_ex = MacroTools.prewalk(s->_replace_me(mc, s), ex)
    eval(Dolo, :(
    let
        $new_ex
    end))
end

eval_with(mc::ModelCalibration, s::AbstractString) = eval_with(mc, _to_expr(s))
eval_with(mc::ModelCalibration, s::Symbol) = _replace_me(mc, s)
eval_with(mc::ModelCalibration, x::Number) = x
eval_with(mc::ModelCalibration, x::AbstractArray) = map(y->eval_with(mc, y), x)

function eval_with(mc::ModelCalibration, d::Associative)
    out = Dict{Symbol,Any}()
    for (k, v) in d
        sk = Symbol(k)
        if sk == :tag
            out[sk] = v
        else
            out[sk] = eval_with(mc, v)
        end
    end
    out
end
