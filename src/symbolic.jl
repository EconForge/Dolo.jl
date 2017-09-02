# -------------- #
# Symbolic Model #
# -------------- #
#
function _get_args(sm, spec)
    # get args
    args = OrderedDict{Symbol,Vector{Tuple{Symbol,Int}}}()
    for (grp, shift, nm) in spec[:eqs]
        grp == "parameters" && continue
        syms = sm.symbols[Symbol(grp)]
        args[Symbol(nm)] = Tuple{Symbol,Int}[(v, shift) for v in syms]
    end
    args
end

function Dolang.FunctionFactory(sm::Model, func_nm::Symbol)
    spec = RECIPES[:dtcc][:specs][func_nm]
    eqs = sm.equations[func_nm]

    # get targets
    target = get(spec, :target, [nothing])[1]
    has_targets = !(target === nothing)
    if has_targets
        target_date = spec[:target][2]
        targets = [Dolang.normalize((s, target_date)) for s in sm.symbols[Symbol(target)]]
    else
        targets =  Symbol[]
    end
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

"""
    build_definition_function(model::Model, defs::Associative=model.definitions)

Generate the code for a Julia function that will evalaute all the definitions
for the model, given an AxisArray containing all states, controls, and
exogenous values.

The signature of the resulting method will be

    evaluate_definitions(::Model, data::AxisArray, p)

where `p` is the vector of model parameters.

The return value is an AxisArray containing the definitions. Will be padded
with NaNs if needed because of non-zero timing in definitions

"""
function build_definition_function{T<:Model}(
        model::T,
        defs::Associative=model.definitions;
        funname::Symbol=:evaluate_definitions
    )
    defs = OrderedDict(defs)
    # used to unpack parameters
    a_factory = first(model.factories)[2]

    # get incidence table for each definition -- skip tracking parameters
    it = Dolang.IncidenceTable(collect(values(defs)), model.symbols[:parameters])

    # need to construct a block for unpacking variables from AxisArrays
    unpack_data_block = Expr(:block)
    for (var, dates) in it.by_var
        _var = Base.QuoteNode(var)
        if haskey(defs, var)
            continue
        end
        for date in dates
            lhs = Dolang.normalize((var, date))

            # if date is non-zero, need to pad array with nans
            if date < 0
                # negative dates requires padding the front of the array
                rhs = :(vcat(fill(NaN, $(abs(date))), data[:, $_var, 1:(capT- $date)]))
            elseif date > 0
                # positive date pad the end
                rhs = :(vcat(data[:, $_var, (1+$date):capT], fill(NaN, $(date))))
            else
                rhs = :(data[:, $_var, :])
            end
            push!(unpack_data_block.args, :($lhs = $rhs))
        end
    end

    # compute incidence one more time, but don't track stuff we already know
    known_vars = vcat(model.symbols[:parameters], get_variables(model))
    it_defs = Dolang.IncidenceTable(
        collect(values(defs)),
        known_vars
    )

    # figure out the order for computing definitions
    def_order = Dolang.solution_order(defs, it_defs, known_vars)
    dot_call(ex) = Expr(:macrocall, Symbol("@__dot__"), Dolang.normalize(ex))
    def_exprs =[Expr(:(=), Dolang.normalize((lhs, 0)), dot_call(rhs))
        for (lhs, rhs) in defs
    ][def_order]
    def_block = Expr(:block, def_exprs...)

    # finally package up the definitions
    # need to reshape all the defs we just computed so we can cat them in dim 2
    cat_expr = Expr(
        :call, :cat, 2,
        [Expr(:call, :reshape, Dolang.normalize((v, 0)), Expr(:tuple, :capN, 1, :capT))
         for v in keys(defs)]...
    )
    # need an axis for the variable dimension
    def_axis = :(Axis{:V}($(collect(keys(defs)))))

    # shove the output into an AxisArray
    out_block = Expr(:block,
        :(out_raw = $cat_expr),
        :(return AxisArray(out_raw, data.axes[1], $(def_axis), (data.axes[3])))
    )

    quote
        function $(funname){T1}(::$T, data::AxisArray{T1,3}, p)
            capT = size(data, 3)
            capN = size(data, 1)
            $(Dolang.param_block(a_factory))
            $unpack_data_block
            $def_block
            $out_block

        end
        function $(funname){T1}(model::$T, tab::AxisArray{T1,2}, p)
            # add leading dimension for N
            data = AxisArray(reshape(tab.data, (1, size(tab)...)), Axis{:N}(1:1), tab.axes...)

            # call method above that does panel
            def_data = $(funname)(model, data, p)

            # drop first dimension from the panel
            return def_data[1, :, :]
        end
    end
end
