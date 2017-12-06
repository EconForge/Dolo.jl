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

############################
# norms for list of points #
############################
import Base.maxabs

function Base.maxabs(s::ListOfPoints)
    t = 0.0
    for p in s
        t = max(t, maximum(p))
    end
    t
end

Base.maxabs(x::Vector{<:ListOfPoints}) = maximum(maxabs.(x))

############################
# serial multiplications #
############################

function invert!(A::Vector)
    # A[i] <- (A[i])^(-1)
    N = length(A)
    for n=1:N
        A[n] = inv(A[n])
    end
end

function premult!(A,B)
    # B[i] <- A[i]*B[i]
    N = length(A)
    for n=1:N
        B[n] = A[n]*B[n]
    end
end

function addmul!(O,A,B)
    # O[i] <- A[i]*B[i]
    N = length(A)
    for n=1:N
        O[n] += A[n]*B[n]
    end
end

##############################
# Conversion from old format #
##############################

vector_to_matrix(v) = Matrix(vec(v)')
# vector_to_matrix(v::Vector) = Matrix(v')
# vector_to_matrix(v::RowVector) = Matrix(v)

#(compat)
to_LOP(::Type{Point{d}}, mat::Matrix) where d = reinterpret(Point{d}, mat', (size(mat,1),))
to_LOP(::Type{Point{d}}, mat::AbstractMatrix) where d = to_LOP(Point{d}, Array(mat))
to_LOP(mat::AbstractArray) = to_LOP(Point{size(mat,2)} ,mat)

function from_LOP(lop)
    d = length(lop[1])
    N = length(lop)
    return reinterpret(Float64, lop, (d,N))'
end

function to_LOJ(mat::Array{Float64,3})
    # list of jacobians
    N,d = size(mat)
    reinterpret(SMatrix{d,d,Float64,d*d}, permutedims(mat,[2,3,1]), (N,))
end
