# ----------------- #
# Parsing utilities #
# ----------------- #
call_expr(var, n) = n == 0 ? symbol(var) :
                             symbol(string(var, "_", n > 0 ? "_" : "m", abs(n)))

function eq_expr(ex::Expr, targets::Union{Vector{Expr},Vector{Symbol}}=Symbol[])
    if isempty(targets)
        return Expr(:call, :(.-), _parse(ex.args[2]), _parse(ex.args[1]))
    end

    # ensure lhs is in targets
    if !(ex.args[1] in targets)
        msg = string("Expected expression of the form `lhs = rhs` ",
                     "where `lhs` is one of $(targets)")
        error(msg)
    end

    Expr(:(=), _parse(ex.args[1]), _parse(ex.args[2]))
end

_parse(x::Symbol) = symbol(string(x, "_"))
_parse(x::Number) = x

function _parse(ex::Expr; targets::Union{Vector{Expr},Vector{Symbol}}=Symbol[])
    @match ex begin
        # translate lhs = rhs  to rhs - lhs
        $(Expr(:(=), :__)) => eq_expr(ex, targets)

        # translate x(n) --> x__n_ and x(-n) -> x_mn_
        var_(shift_Integer) => _parse(call_expr(var, shift))

        # + and * can take multiple arguments. Need to slurp them
        +(a_, b_) => Expr(:call, :(.+), _parse(a), _parse(b))
        +(a_, b_, c__) => _parse(Expr(:call, :(+), Expr(:call, :(.+), a, b), c...))
        *(a_, b_) => Expr(:call, :(.*), _parse(a), _parse(b))
        *(a_, b_, c__) => _parse(Expr(:call, :(*), Expr(:call, :(.*), a, b), c...))

        # -, /, ^ can only take two at a time
        -(a_, b_) => Expr(:call, :(.-), _parse(a), _parse(b))
        /(a_, b_) => Expr(:call, :(./), _parse(a), _parse(b))
        ^(a_, b_) => Expr(:call, :(.^), _parse(a), _parse(b))

        # Other func calls. Just parse args. Allows arbitrary Julia functions
        f_(a__) => Expr(:call, f, map(_parse, a)...)

        # the bottom, just insert numbers and symbols
        x_Symbol_Number => x

        _ => error("Not sure what I just saw")
    end
end

_parse(s::AbstractString; kwargs...) = _parse(parse(s); kwargs...)

# -------------------- #
# Handling Definitions #
# -------------------- #

"Returns `true` if `ex` has the form `x(i::Int)` -- false otherwise"
_is_time_shift(ex::Expr) = ex.head == :call &&
                           length(ex.args) == 2 &&
                           isa(ex.args[2], Int)

"""
Returns `true` if `ex` has the form `x(i::Int)` and `x` is a variable in `sm`
-- false otherwise
"""
_is_time_shift(sm::ASM, ex::Expr) = _is_time_shift(ex) &&
                                    haskey(sm.calibration, ex.args[1])

"""
in `ex` replace all occurances of `sym(i)` with `sym(i+shift)`

If `sym` is found on its own, replace it with `sym(shift)`
"""
function _replace_with_shift(sm::ASM, ex::Expr, sym::Symbol, shift::Int)
    # if this is a time shift itself, apply additional shift
    if _is_time_shift(sm, ex)
        sym, pre_shift = eq.args
        return Expr(:call, sym, pre_shift+shift)
    end

    # otherwise just go through the rest of the args
    Expr(ex.head, [_replace_with_shift(sm, _, sym, shift) for _ in ex.args]...)
end

"Numbers just go through"
_replace_with_shift(sm::ASM, n::Number, sym::Symbol, shift::Int) = n

"If the symbol is what we're looking for, add shift, otherwise pass through"
_replace_with_shift(sm::ASM, input::Symbol, sym::Symbol, shift::Int) =
    input == sym ? Expr(:call, sym, shift) : input

"""
Replace an occurance of `def(i)` with the body of the definition `def` at time
shift `i`.
"""
_call_definition(sm::ASM, defn::Expr, shift::Int=0) =
    Expr(defn.head, [_call_definition(sm, _, shift) for _ in defn.args]...)

"""
`_call_definition(sm::ASM, sym::Symbol, shift::Int=0)`

Do one of 3 things:

1. If `sym` is in `sm.definitions`, then apply the shift to the contents of
`sm.definitions[sym]`.
2. if `sym` is somewhere in `sm.symbols` (excluding sm.symbol[:sparameters]) apply
the time sift.
3. If it is any other symbol, just let it through
"""
function _call_definition(sm::ASM, sym::Symbol, shift::Int=0)
    if haskey(sm.definitions, sym)
        return _call_definition(sm, sm.definitions[sym], shift)
    end

    for grp in keys(sm.symbols)
        grp == :parameters && continue
        for v in sm.symbols[grp]
            sym == v && return :($(v)($shift))
        end
    end

    # otherwise just return the symbol
    return sym
end

"Let numbers through"
_call_definition(sm::ASM, x::Number, shift::Int) = x

"the bottom, just let everything else through"
_apply_definitions(sm::ASM, x) = x

"If `sym` is a definition, recursively resolve, otherwise let `sym` through"
_apply_definitions(sm::ASM, sym::Symbol) =
    haskey(sm.definitions, sym) ? _call_definition(sm, sm.definitions[sym]): sym

"""
`_apply_definitions(sm::ASM, eq::Expr)`

Do one of two things:

1. if `eq` is of the form `def(i::Int)` and `def` is a key in `sm.definitions`:
apply a time shifted version of `sm.definitions[def]`
2. Otherwise, construct the exact same expression, recursively calling this
function for all items in `eq.args`
"""
function _apply_definitions(sm::ASM, eq::Expr)
    if _is_time_shift(eq)
        sym = eq.args[1]
        shift = eq.args[2]
        if haskey(sm.definitions, sym)
            defn = sm.definitions[sym]
            return _call_definition(sm, defn, shift)
        end
    end

    # in all other cases (including other :call exprs) just apply to all args
    Expr(eq.head, [_apply_definitions(sm, _) for _ in eq.args]...)
end

"Map `_apply_definitions` to each element in the vector"
_apply_definitions(sm::ASM, eqs::Vector{Expr}) =
    map(_ -> _apply_definitions(sm, _), eqs)

# -------- #
# Compiler #
# -------- #
@inline _unpack_var(x::AbstractVector, i::Integer) = x[i]
@inline _unpack_var(x::AbstractMatrix, i::Integer) = sub(x, :, i)

function _param_block(sm::ASM)
    params = sm.symbols[:parameters]
    Expr(:block,
         [:($(_parse(params[i])) = _unpack_var(p, $i)) for i in 1:length(params)]...)
end

function _single_arg_block(sm::ASM, arg_name::Symbol, arg_type::Symbol,
                           shift::Int, Ndim::Int=1)
    if arg_name == :p && arg_type == :parameters
        @assert shift == 0
        return _param_block(sm)
    end

    nms = sm.symbols[arg_type]
    # TODO: extract columns at a time when Ndim > 1
    Expr(:block,
         [:($(_parse("$(nms[i])($(shift))")) = _unpack_var($(arg_name),$i))
            for i in 1:length(nms)]...)
end

@inline _assign_var(lhs::AbstractVector, rhs::Number, i) = setindex!(lhs, rhs, i)
@inline _assign_var(lhs::AbstractMatrix, rhs::AbstractVector, i) = setindex!(lhs, rhs, :, i)

"returns an expression `:(lhs[i] = rhs)`"
_assign_var_expr(lhs, rhs, i) = :(_assign_var($lhs, $rhs, $i))

"Evaluates main expressions in a function group and fills `out` with results"
function _main_body_block(sm::ASM, targets::Vector{Symbol}, exprs::Vector{Expr})
    n_expr = length(exprs)
    parsed_exprs = map(x->_parse(x; targets=targets), exprs)
    if isempty(targets)
        # no targets, just set out = parsed_expr
        assignments = map((rhs,i)->_assign_var_expr(:out, rhs, i),
                          parsed_exprs, 1:n_expr)
        func_block = Expr(:block, assignments...)
    else
        # otherwise, need to parse the targets, evaluate them, and then set
        # elements of out equal to the targets
        parsed_targets = map(_parse, targets)
        assignments = map((rhs,i)->_assign_var_expr(:out, rhs, i),
                          parsed_targets, 1:n_expr)
        func_block = Expr(:block, parsed_exprs..., assignments...)
    end
end

"checks that lhs of expr is equal to target"
function _check_target(expr::Expr, target::Symbol, func_nm::Symbol, i::Int)
    out = @match expr begin
        $(Expr(:(=), :__)) => expr.args[1]
        _ => false
    end

    if isa(out, Bool)
        error("Equations of type $func_nm must be of the form `lhs=rhs`")
    end

    # otherwise we got a symbol
    if out != target
        msg = string("In equation of type $func_nm expected $target on the ",
                     " left hand side of the equation $(i)")
        error(msg)
    end

end

function _check_targets(exprs, targets, func_nm)
    nt = length(targets)
    ne = length(exprs)
    if nt != ne
        msg = "For equation group $func_nm expected $nt equations, found $ne"
        error(msg)
    end

    map((e, t, i) -> _check_target(e, t, func_nm, i), exprs, targets, 1:nt)
    nothing
end

_output_size(n_expr::Int, args::AbstractVector...) = (n_expr,)

function _output_size(n_expr::Int, args...)
    n_row = 0
    # get maximum number of rows among matrix arguments
    for a in args
        if isa(a, AbstractMatrix)
            nr = size(a, 1)
            if n_row > 0
                # we already found another matrix argument, now we can enforce
                # that all matrix arguments have conformable shapes
                if nr != n_row
                    msg = string("Unconformable argument sizes. For vectorized ",
                                 "evaluation all matrix arguments must have ",
                                 "the same number of rows.")
                    throw(DimensionMismatch(msg))
                end
            else
                # we need to update n_row
                n_row = nr
            end
        end
    end
    (n_row, n_expr)
end

_allocate_out(T::Type, n_expr::Int, args::AbstractVector...) = Array(T, n_expr)

function _allocate_out(T::Type, n_expr::Int, args...)
    sz = _output_size(n_expr, args...)
    Array(T, sz[1], sz[2])
end

function compile_equation(sm::ASM, func_nm::Symbol; print_code::Bool=false)
    # extract spec from recipe
    spec = RECIPES[model_type(sm)][:specs][func_nm]

    # get expressions from symbolic model
    exprs = sm.equations[func_nm]

    numeric_mod = _numeric_mod_type(sm)

    bang_func_nm = symbol(string(func_nm), "!")

    if length(exprs) == 0
        msg = "Model did not specify functions of type $(func_nm)"
        code = quote
            function $(func_nm)(::$(numeric_mod), args...)
                error($msg)
            end

            function $(bang_func_nm)(::$(numeric_mod), args...)
                error($msg)
            end
        end
        return code
    end

    # first apply definitions
    exprs = _apply_definitions(sm, exprs)

    # extract information from spec
    target = get(spec, :target, [nothing])[1]
    eqs = spec[:eqs]  # required, so we don't provide a default
    non_aux = filter(x->x[1] != "auxiliaries", eqs)
    arg_names = Symbol[symbol(x[3]) for x in non_aux]
    arg_types = Symbol[symbol(x[1]) for x in non_aux]
    arg_shifts = Int[x[2] for x in non_aux]
    only_aux = filter(x->x[1] == "auxiliaries", eqs)
    aux_shifts = Int[x[2] for x in only_aux]

    # extract targets and make sure they appear in the correct order in exprs
    has_targets = !(target === nothing)
    targets = has_targets ? sm.symbols[symbol(target)] : Symbol[]
    has_targets && _check_targets(exprs, targets, func_nm)

    # build function block by block
    all_arg_blocks = map((a,b,c) -> _single_arg_block(sm, a, b, c),
                         arg_names, arg_types, arg_shifts)
    arg_block = Expr(:block, vcat([b.args for b in all_arg_blocks]...)...)
    all_aux_blocks = map(n -> _aux_block(sm, n), aux_shifts)
    aux_block = Expr(:block, vcat([b.args for b in all_aux_blocks]...)...)
    main_block = _main_body_block(sm, targets, exprs)

    # construct the body of the function
    body = quote
        $(arg_block)
        $(aux_block)
        $(main_block)
        out  # return out
    end

    typed_args = [Expr(:(::), s, :(Union{AbstractVector,AbstractMatrix}))
                  for s in arg_names]

    # build the new type and implement methods on Base.call that we need
    code = quote
        # non-allocating function
        function $(bang_func_nm)(out, ::$(numeric_mod), $(typed_args...))
            expected_size = _output_size($(length(exprs)), $(arg_names...))
            if size(out) != expected_size
                msg = "Expected out to be size $(expected_size), found $(size(out))"
                throw(DimensionMismatch(msg))
            end
            $body  # evaluates equations and populates `out`
        end

        # allocating version
        function $(func_nm)(mod::$(numeric_mod), $(typed_args...))
            out = _allocate_out(eltype($(arg_names[1])),
                                $(length(exprs)), $(arg_names...))
            $(bang_func_nm)(out, mod, $(arg_names...))
        end

        # TODO: can we use broadcast! to get pretty far towards guvectorize?
    end

    print_code && println(code)

    code

end
