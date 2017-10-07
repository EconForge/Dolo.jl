struct BasisMatricesDR{S,T,nx,Titp<:Vector{<:BM.Interpoland}} <: AbstractDecisionRule{S,T,nx}
    grid_exo::S
    grid_endo::T
    itp::Titp
end

nodes(dr::BasisMatricesDR) = BM.nodes(dr.itp[1].basis)[1]
Base.ndims(dr::BasisMatricesDR) = ndims(dr.itp[1].basis)
n_nodes(dr::BasisMatricesDR) = length(dr.itp[1].basis)

# grid_exo <: EmptyGrid
function BasisMatricesDR(
        grid_exo::S, grid_cart::T, ::Union{Val{nx},Type{Val{nx}}},
        param_type::Type{<:BM.BasisParams}=BM.ChebParams,
        args::Tuple=tuple()
    ) where S <: EmptyGrid where T <: CartesianGrid where nx
    params = (param_type.(grid_cart.n, grid_cart.min, grid_cart.max, args...)...)
    basis = BM.Basis(params)
    coefs = [zeros(length(basis), nx)]
    itp = [BM.Interpoland(basis, c) for c in coefs]
    BasisMatricesDR{S,T,nx,typeof(itp)}(grid_exo, grid_cart, itp)
end

# grid_exo <: UnstructuredGrid
function BasisMatricesDR(
        grid_exo::S, grid_cart::T,
        ::Union{Val{nx},Type{Val{nx}}},
        param_type::Type{<:BM.BasisParams}=BM.ChebParams,
        args::Tuple=tuple()
    ) where S <: UnstructuredGrid where T <: CartesianGrid where nx
    params = (param_type.(grid_cart.n, grid_cart.min, grid_cart.max, args...)...)
    basis = BM.Basis(params)
    coefs = [zeros(length(basis), nx) for i in 1:n_nodes(grid_exo)]
    itp = [BM.Interpoland(basis, c) for c in coefs]
    BasisMatricesDR{S,T,nx,typeof(itp)}(grid_exo, grid_cart, itp)
end

## grid_exo <: Union{EmptyGrid,UnstructuredGrid}
function set_values!(
        dr::BasisMatricesDR{<:G}, values::Vector{Matrix{Float64}}
    ) where G <: Union{EmptyGrid,UnstructuredGrid}
    for i in 1:length(values)
        BM.update_coefs!(dr.itp[i], values[i])
    end
end

function set_values!(
        dr::BasisMatricesDR{<:G,<:Grid,nx},
        values::Vector{<:Array{Value{nx}}}
    ) where G <: Union{EmptyGrid,UnstructuredGrid} where nx
    if length(values) != length(dr.itp)
        msg = "The length of values ($(length(values))) is not the same "
        msg *= "as the length of the coefficient Vector ($(length(dr.itp)))"
        error(msg)
    end

    for i in 1:length(values)
        N = length(values[i])
        data = reinterpret(Float64, values[i], (nx, N))'
        BM.update_coefs!(dr.itp[i], data)
    end
end

## grid_exo <: EmptyGrid
evaluate(dr::BasisMatricesDR{EmptyGrid}, z::Matrix) = dr.itp[1](z)

## grid_exo <: UnstructuredGrid
function evaluate(dr::BasisMatricesDR{<:UnstructuredGrid}, i::Int, z::AbstractMatrix)
    @boundscheck begin
        n_funcs = length(dr.itp)
        if i > n_funcs
            msg = "Only $n_funcs are known, but function $i was requested"
            throw(BoundsError(msg))
        end
    end
    dr.itp[i](z)
end

## fully general
function evaluate(dr::BasisMatricesDR{<:Grid,<:Grid{d}}, points::Vector{Point{d}}) where d
    N = length(points)
    mat = reinterpret(Float64, points, (d, N))'
    evaluate(dr, mat)
end

function evaluate(dr::BasisMatricesDR{<:Grid,<:Grid{d}}, i::Int, points::Vector{Point{d}}) where d
    N = length(points)
    mat = reinterpret(Float64, points, (d, N))'
    evaluate(dr, i, mat)
end
