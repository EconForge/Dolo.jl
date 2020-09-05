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


# ----- #
# Tools #
# ----- #


# Used for extracting complementarity conditions from arbitrage equations
function _handle_arbitrage(arb, controls)

    controls_lb = Expression[]
    controls_ub = Expression[]
    arbitrage = Expression[]
    for (i, v) in enumerate(arb)

        eq = Dolang.parse_string(v)

        parts = Dolang.match_equation(eq)
        if length(parts) == 1
            push!(arbitrage, parts[1])
            push!(controls_lb, :(-Inf))
            push!(controls_ub, :Inf)
        elseif length(parts) == 2
            f = :($(parts[2])-$(parts[1]))
            push!(arbitrage, f)
            push!(controls_lb, :(-Inf))
            push!(controls_ub, :Inf)
        else
            # length(parts) == 5
            f = parts[1]

            push!(arbitrage, f)
            push!(controls_lb, parts[2])
            push!(controls_ub, parts[4])

            # verify that the control is what we want
            mid = parts[3]
            if mid != controls[i]
                msg = string("Error in complementarity condition. ",
                                "Expected $(controls[i]) found $mid")
                error(msg)
            end
        end

    end
    controls_lb, controls_ub, arbitrage
end

