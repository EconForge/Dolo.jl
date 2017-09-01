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


function to_LOP(mat::Matrix{Float64})
    N,d = size(mat)
    reinterpret(Point{d}, mat', (N,))
end

function to_LOJ(mat::Array{Float64})
    # list of jacobians
    N,d = size(mat)
    reinterpret(SMatrix{d,d,Float64,d*d}, permutedims(mat,[2,3,1]), (N,))
end


function from_LOP(lop)
    d = length(lop[1])
    N = length(lop)
    return reinterpret(Float64, lop, (d,N))'
end
