_expr_or_number(x::Union{AbstractString,Symbol,Expr}) = _to_expr(x)
_expr_or_number(x::Number) = x

_to_Float64(x::Real) = convert(Float64, x)
_to_Float64(x::AbstractArray) = map(Float64, x)

appenddim(A) = reshape(A, size(A)..., 1)

# does kronecker product of matrices according to fortran order
function fkron(A::AbstractMatrix, B::AbstractMatrix)
    kron(B, A)
end

rmerge(default_struct, add_options) = add_options
function rmerge(def_s::Associative, add_s::Associative)
    kl = collect(keys(def_s))
    kr = collect(keys(add_s))
    resp = Dict()
    for k in intersect(kl,kr)
        resp[k] = rmerge(def_s[k],add_s[k])
    end
    for k in setdiff(kl, kr)
        resp[k] = def_s[k]
    end
    for k in setdiff(kr, kl)
        resp[k] = add_s[k]
    end
    return resp
end
