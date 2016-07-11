_replace_star_star(s::AbstractString) = replace(s, "**", "^")

_to_expr(x::Expr) = x
_to_expr(x::Union{Symbol,Number}) = Expr(:block, x)
_to_expr(x::AbstractString) = _to_expr(parse(_replace_star_star(x)))

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

_to_Float64(x::Real) = convert(Float64, x)
_to_Float64(x::AbstractArray) = map(Float64, x)

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
                    sol = eval(Dolo,
                               :(let
                                    $([:($x=$y) for (x, y) in solutions]...);
                                    $expr
                                 end)
                               )
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
