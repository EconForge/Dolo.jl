# ----------------- #
# Parsing utilities #
# ----------------- #
call_expr(var, n) = n == 0 ? symbol(var) :
                             symbol(string(var, "_", n > 0 ? "_" : "m", abs(n)))

function eq_expr(ex::Expr, targets::Union{Vector{Expr},Vector{Symbol}}=Symbol[])
    if isempty(targets)
        return Expr(:call, :(-), _parse(ex.args[2]), _parse(ex.args[1]))
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

        # Other func calls. Just parse args. Allows arbitrary Julia functions
        f_(a__) => Expr(:call, f, map(_parse, a)...)

        # the bottom, just insert numbers and symbols
        x_Symbol_Number => x

        _ => error("Not sure what I just saw")
    end
end

_parse(s::AbstractString; kwargs...) = _parse(parse(s); kwargs...)

# -------- #
# Compiler #
# -------- #

function _param_block(sm::ASM)
    params = sm.symbols[:parameters]
    Expr(:block,
         [:(@inbounds $(_parse(params[i])) = p[$i]) for i in 1:length(params)]...)
end

function _aux_block(sm::ASM, shift::Int)
    target = RECIPES[model_spec(sm)][:specs][:auxiliary][:target][1]
    targets = sm.symbols[symbol(target)]
    exprs = sm.equations[:auxiliary]

    if shift == 0
        return Expr(:block, [_parse(ex; targets=targets) for ex in exprs]...)
    end

    # aux are strict functions of states and controls at time 1. We can do some
    # string manipulation of each expr and replace each with the shifted
    # version
    string_expr = map(string, exprs)

    # now we have to deal with shift
    for grp in (:states, :controls, :auxiliaries)
        for i in 1:length(string_expr)
            for sym in sm.symbols[grp]
                pat = Regex("\\b$(sym)\\b")
                rep = "$sym($shift)"
                string_expr[i] = replace(string_expr[i], pat, rep)
            end
        end
    end

    # now adjust targets. Don't need anything fancy b/c we know we have a
    # symbol
    targets = [:($(t)($(shift))) for t in targets]

    Expr(:block, [_parse(ex; targets=targets) for ex in string_expr]...)
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
         [:(@inbounds $(_parse("$(nms[i])($(shift))")) = $(arg_name)[$i])
            for i in 1:length(nms)]...)
end

"returns an expression `:(lhs[i] = rhs)`"
_assign_single_el(lhs, rhs, i) = :($lhs[$i] = $rhs)

"Evaluates main expressions in a function group and fills `out` with results"
function _main_body_block(sm::ASM, targets::Vector{Symbol}, exprs::Vector{Expr})
    n_expr = length(exprs)
    parsed_eprs = map(x->_parse(x; targets=targets), exprs)
    assignments = map((rhs,i)->_assign_single_el(:out, rhs, i),
                      parsed_eprs, 1:n_expr)
    func_block = Expr(:block, assignments...)
end

function compile_equation(sm::ASM, func_nm::Symbol)
    # extract spec from recipe
    spec = RECIPES[model_spec(sm)][:specs][func_nm]

    # generate a new type name
    tnm = gensym(func_nm)

    # get expressions from symbolic model
    exprs = sm.equations[func_nm]

    if length(exprs) == 0
        # we are not able to use this equation type. Just create a dummy type
        # and function that throws an error explaining what went wrong
        msg = "Model did not specify functions of type $(func_nm)"
        code = quote
            immutable $tnm <: AbstractDoloFunctor
            end
            function evaluate(::$tnm, args...)
                error($msg)
            end

            function evaluate!(::$tnm, args...)
                error($msg)
            end

            $tnm()  # see note below
        end
        return code
    end

    # extract information from spec
    target = get(RECIPES[model_spec(sm)][:specs][func_nm], :target, [nothing])[1]
    targets = target === nothing ? Symbol[] : sm.symbols[symbol(target)]
    eqs = spec[:eqs]  # required, so we don't provide a default
    non_aux = filter(x->x[1] != "auxiliaries", eqs)
    arg_names = Symbol[symbol(x[3]) for x in non_aux]
    arg_types = Symbol[symbol(x[1]) for x in non_aux]
    arg_shifts = Int[x[2] for x in non_aux]
    only_aux = filter(x->x[1] == "auxiliaries", eqs)
    aux_shifts = Int[x[2] for x in only_aux]

    # build function block by block
    all_arg_blocks = map((a,b,c) -> _single_arg_block(sm, a, b, c),
                         arg_names, arg_types, arg_shifts)
    arg_block = Expr(:block, all_arg_blocks...)
    all_aux_blocks = map(n -> _aux_block(sm, n), aux_shifts)
    aux_block = Expr(:block, all_aux_blocks...)
    main_block = _main_body_block(sm, targets, exprs)

    # construct the body of the function
    body = quote
        $(arg_block)
        $(aux_block)
        $(main_block)
        out  # return out
    end

    typed_args = [Expr(:(::), s, :AbstractVector) for s in arg_names]

    # build the new type and implement methods on Base.call that we need
    code = quote
        immutable $tnm <: AbstractDoloFunctor
        end

        # non-allocating function
        function evaluate!(::$tnm, $(typed_args...), out)
            $body  # evaluateuates equations and populates `out`
        end

        # allocating version
        function evaluate(o::$tnm, $(typed_args...))
            out = Array(eltype($(arg_names[1])), $(length(exprs)))
            evaluate!(o, $(arg_names...), out)
        end

        # last line of this block is the singleton instance of the type
        # This means you should do `obj = eval(code)`
        $tnm()
        # TODO: can we use broadcast! to get pretty far towards guvectorize?
    end

    code

end
