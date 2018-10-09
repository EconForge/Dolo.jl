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
function rmerge(def_s::AbstractDict, add_s::AbstractDict)
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
# import Base.maxabs
#
# function maxabs(s::ListOfPoints)
#     t = 0.0
#     for p in s
#         t = max(t, maximum(p))
#     end
#     t
# end


function maxabs(s::AbstractVector{Point{d}}) where d
    t = 0.0
    for p in s
        t = max(t, maximum(p))
    end
    t
end

maxabs(x::AbstractVector{<:AbstractVector{Point{d}}}) where d= maximum(maxabs.(x))

############################
# serial multiplications #
############################

function invert!(A::AbstractVector)
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
to_LOP(::Type{Point{d}}, mat::Matrix) where d = reshape(reinterpret(Point{d}, vec(copy(mat'))), (size(mat,1),))
to_LOP(::Type{Point{d}}, mat::AbstractMatrix) where d = to_LOP(Point{d}, Array(mat))
to_LOP(mat::AbstractArray) = to_LOP(Point{size(mat,2)} ,mat)

function from_LOP(lop)
    d = length(lop[1])
    N = length(lop)
    return reshape(reinterpret(Float64, vec(lop)), (d,N))'
end

function to_LOJ(mat::Array{Float64,3})
    # list of jacobians
    N,d = size(mat)
    reshape(reinterpret(SMatrix{d,d,Float64,d*d}, vec(permutedims(mat,[2,3,1]))), (N,))
end
