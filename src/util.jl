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

function solution_order(d::OrderedDict, it::Dolang.IncidenceTable)
    # unpack some data
    vars = collect(keys(d))
    n = length(d)

    # allocate
    out = zeros(Int, length(d))

    # Start with indices for equations that are purely numerical
    front = setdiff(1:n, keys(it.by_eq))
    out[front] = 1:length(front)
    solved = vars[front]
    to_solve = deepcopy(it.by_eq)

    # now start stepping through equations
    ix = length(front)
    for _ in 2:n
        for (eq, eq_vars) in to_solve
            can_solve = true
            for (var, dates) in eq_vars
                if !in(var, solved)
                    can_solve = false
                    break
                end
            end
            if can_solve
                out[eq] = ix+=1
                push!(solved, vars[eq])
                pop!(to_solve, eq)
            end
        end
    end

    !isempty(to_solve) && error("Not triangular system")

    return sortperm(out)
end


function solution_order(_d::Associative)
    d = OrderedDict(_d)
    it = Dolang.IncidenceTable(collect(values(d)))
    solution_order(d, it)
end

solve_triangular_system(sm::ASM) = solve_triangular_system(sm.calibration)
solve_triangular_system(d::Associative) = solve_triangular_system(OrderedDict(d))

function solve_triangular_system(d::OrderedDict)
    sol_order = solution_order(d)

    # extract expressions and variable names in proper order
    nms = collect(keys(d))[sol_order]
    exprs = collect(values(d))[sol_order]

    # build expression to evaluate system in correct order
    to_eval = Expr(:block)
    to_eval.args = [:($(i[1])=$(i[2])) for i in zip(nms, exprs)]

    # add one line to return a tuple of all data
    ret = Expr(:tuple); ret.args = nms

    # now evaluate and get data
    data = eval(Dolo, :(let
                        $to_eval;
                        $ret
                        end))

    OrderedDict(zip(nms, data))
end
