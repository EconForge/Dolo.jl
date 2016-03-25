_to_expr(x::Expr) = x
_to_expr(x::Union{Symbol,Number}) = Expr(:block, x)
_to_expr(x::AbstractString) = _to_expr(parse(x))

_expr_or_number(x::Union{AbstractString,Symbol,Expr}) = _to_expr(x)
_expr_or_number(x::Number) = x

inf_to_Inf(x::Number) = x
inf_to_Inf(x::Symbol) = @match x begin
    inf => Inf
    _ => x
end

inf_to_Inf(x::Expr) = @match x begin
    -inf => -Inf
    -x_Number => -x
    f_(a__) => Expr(:call, f, map(inf_to_Inf, a)...)
end

solve_triangular_system(sm::ASM) = solve_triangular_system(sm.calibration)

function solve_triangular_system(dict::Associative)
    solutions = Dict{Symbol,Float64}()
    finished = false
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

    length(solutions) < length(dict) &&  error("Not a triangular system")

    # reorder solutions to match sm.calibration
    OrderedDict{Symbol,Float64}([(k, solutions[k]) for k in keys(dict)])
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
