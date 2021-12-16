function label_density(μ, symbols,  grid_exo::EmptyGrid, grid_endo::UCGrid)
    endo_names = symbols[:states]
    sc = scales(grid_endo)
    d = Dict(endo_names[i]=>sc[i] for i=1:length(grid_endo.n))
    return AxisArray(μ; d...)
end


function ergodic_distribution(model, sol)
    return ergodic_distribution(model, sol.dr, sol.dr.grid_exo, sol.dr.grid_endo, sol.dprocess)
end

function ergodic_distribution(model::Model{T, Q}, sol) where T where Q<:Dolo.IIDExogenous
    P, μ = ergodic_distribution(model, sol.dr, sol.dr.grid_exo, sol.dr.grid_endo, sol.dprocess)
    
    return label_density(μ, model.symbols,  sol.dr.grid_exo, sol.dr.grid_endo)
end


function ergodic_distribution(model, dr, exo_grid:: UnstructuredGrid, endo_grid:: UCGrid, dp)
    N_m = n_nodes(exo_grid)
    N_s = n_nodes(endo_grid)
    N = N_m*N_s
    Π = zeros(N_m, N_s, N_m, endo_grid.n...)
    s = nodes(endo_grid)
    p = SVector(model.calibration[:parameters]...)
    a = SVector(endo_grid.min...)
    b = SVector(endo_grid.max...)
    for i_m in 1:n_nodes(exo_grid)
        x = dr(i_m, s)
        m = node(exo_grid, i_m)
        for i_M in 1:n_inodes(dp, i_m)
            M = inode(Point, dp, i_m, i_M)
            w = iweight(dp, i_m, i_M)
            S = transition(model, m, s, x, M, p)
            S = [(S[n]-a)./(b-a) for n=1:length(S)]
            if typeof(Π) <: AbstractArray{Float64,4}
                trembling_hand!(view(Π,i_m, :, i_M, :), S, w)
            elseif typeof(Π) <: AbstractArray{Float64,5}
                trembling_hand!(view(Π,i_m, :, i_M, :, :), S, w)
            end
        end
    end
    Π0 = (reshape(Π,N,N))' - Matrix(I, N, N)
    Π0[end,:] .= 1
    B = zeros(N)
    B[end] = 1.0
    μ = Π0\B

    μ = reshape(μ, N_m, N_s)
    return Π, μ
end

function ergodic_distribution(model, dr, exo_grid:: UCGrid, endo_grid:: UCGrid, dp)
    N_m = n_nodes(exo_grid)
    N_s = n_nodes(endo_grid)
    N = N_m*N_s
    Π = zeros(N_m, N_s, exo_grid.n..., endo_grid.n...)
    s = nodes(endo_grid)
    p = SVector(model.calibration[:parameters]...)
    a = SVector(exo_grid.min..., endo_grid.min...)
    b = SVector(exo_grid.max..., endo_grid.max...)
    for i_m in 1:n_nodes(exo_grid)
        x = dr(i_m, s)
        m = node(exo_grid, i_m)
        for i_M in 1:n_inodes(dp, i_m)
            M = inode(Point, dp, i_m, i_M)
            w = iweight(dp, i_m, i_M)
            S = transition(model, m, s, x, M, p)
            V = [(SVector(M..., el...)-a)./(b.-a) for el in S]
            # return Π, i_m, V, w
            if typeof(Π) <: AbstractArray{Float64,4}
                trembling_hand!(view(Π,i_m, :, :, :), V, w)
            elseif typeof(Π) <: AbstractArray{Float64,3}
                trembling_hand!(view(Π,i_m, :, :), V, w)
            end
        end
    end

    Π0 = (reshape(Π,N,N))' - Matrix(I, N, N)
    Π0[end,:] .= 1
    B = zeros(N)
    B[end] = 1.0
    μ = Π0\B

    μ = reshape(μ, N_m, N_s)
    return reshape(Π, N, N), μ
end


function ergodic_distribution(model, dr, exo_grid:: EmptyGrid, endo_grid:: UCGrid, dp)
    N_m = 1
    N_s = n_nodes(endo_grid)
    N = N_m*N_s
    Π = zeros(N_s, endo_grid.n...)
    s = nodes(endo_grid)
    p = SVector(model.calibration[:parameters]...)
    a = SVector(endo_grid.min...)
    b = SVector(endo_grid.max...)
    i_m = 1
    x = dr(s)
    m = SVector(model.calibration[:exogenous]...)
    for i_M in 1:n_inodes(dp, i_m)
        M = inode(Point, dp, i_m, i_M)
        w = iweight(dp, i_m, i_M)
        S = transition(model, m, s, x, M, p)
        S = [(S[n]-a)./(b-a) for n=1:length(S)]
        trembling_hand!(Π, S, w)
    end
    Π0 = (reshape(Π,N,N))' - Matrix(I, N, N)
    Π0[end,:] .= 1
    B = zeros(N)
    B[end] = 1.0
    μ = Π0\B

    μ = reshape(μ, endo_grid.n...)
    return reshape(Π, N, N), μ
end


"""
Computes the outer product.

# Argument
* `λn_weight_vector::Vararg{Point{2}}`: tuple of Point{2} to be multiplied by outer product

# Returns
* the outer product
"""
function outer(λn_weight_vector::Vararg{Point{2}})
    return [prod(e) for e in Iterators.product(λn_weight_vector...)]
end


"""
Updates A.

# Arguments
* `A`: the transition matrix that will be updated.
* `x::Vector{Point{d}}` : vector of controls.
* `w::Float64` : vector of weights.

# Modifies
* `A` : the updated transition matrix 
"""
function trembling_hand!(A, x::Vector{Point{d}}, w::Float64) where d
    
    @assert ndims(A) == d+1
    shape_A = size(A)
    grid_dimension = d
    δ =  SVector{d,Float64}(1.0./(shape_A[1+i]-1) for i in 1:d )
    

    for n in 1:shape_A[1]

        xn = x[n]
        xn = min.(max.(xn, 0.0),1.0)
        qn = div.(xn, δ)
        qn = max.(0, qn)
        qn = min.(qn, shape_A[2:d+1].-2)
        λn = (xn./δ.-qn) # ∈[0,1[ by construction
        qn_ = round.(Int,qn) + 1
        
        λn_weight_vector = tuple( (SVector(1-λn[i],λn[i]) for i in 1:d)... )

        indexes_to_be_modified = tuple(n, UnitRange.(qn_,qn_.+1)...)

        # Filling transition matrix
        rhs = outer(λn_weight_vector...)
        A[indexes_to_be_modified...] .+= rhs
        
    end

end

"""
Calculates the new transition matrix for a given model, a given discretized exogenous process, given control values (x0) and given grids (exogenous and endogenous).

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `dprocess::`: Discretized exogenous process.
* `x0::Dolo.MSM{SVector{2, Float64}}`: Initial control values.
* `exo_grid`: Exogenous grid that can be of type either UnstructuredGrid or UCGrid or EmptyGrid (in the three following functions).
* `endo_grid::UCGrid`: Endogenous grid.
* `exo`: nothing or (z0, z1)

# Returns
* `Π0::`: New transition matrix.
"""

function new_transition(model, dp, x0, exo_grid:: UnstructuredGrid, endo_grid:: UCGrid; exo=nothing)

    parms = SVector(model.calibration[:parameters]...)

    N_m = n_nodes(exo_grid)
    N_s = n_nodes(endo_grid)
    N = N_m*N_s
    Π = zeros(N_m, N_s, N_m, endo_grid.n...)
    s = nodes(endo_grid)
    a = SVector(endo_grid.min...)
    b = SVector(endo_grid.max...)
    for i_m in 1:n_nodes(exo_grid)
        x = x0.views[i_m]
        m = node(exo_grid, i_m)
        if !(exo === nothing)
            m = Dolo.repsvec(exo[1], m)   # z0
        end
        for i_M in 1:n_inodes(dp, i_m)
            M = inode(Point, dp, i_m, i_M)
            if !(exo === nothing)
                M = Dolo.repsvec(exo[2], M)   # z1
            end
            w = iweight(dp, i_m, i_M)
            S = transition(model, m, s, x, M, parms)
            S = [(S[n]-a)./(b-a) for n=1:length(S)]
            trembling_hand!(view(Π,tuple(i_m,:,i_M,(Colon() for k in 1:(ndims(Π)-3))...)...), S, w)
        end
    end
    Π0 = (reshape(Π,N,N))

    return Π0 
end

function new_transition(model, dp, x0, exo_grid:: UCGrid, endo_grid:: UCGrid; exo=nothing)

    parms = SVector(model.calibration[:parameters]...)

    N_m = n_nodes(exo_grid)
    N_s = n_nodes(endo_grid)
    N = N_m*N_s
    Π = zeros(N_m, N_s, exo_grid.n..., endo_grid.n...)
    s = nodes(endo_grid)
    a = SVector(exo_grid.min..., endo_grid.min...)
    b = SVector(exo_grid.max..., endo_grid.max...)
    for i_m in 1:n_nodes(exo_grid)
        x = x0.views[i_m]
        m = node(exo_grid, i_m)
        if !(exo === nothing)
            m = Dolo.repsvec(exo[1], m)   # z0
        end
        for i_M in 1:n_inodes(dp, i_m)
            M = inode(Point, dp, i_m, i_M)
            if !(exo === nothing)
                M = Dolo.repsvec(exo[2], M)   # z1
            end
            w = iweight(dp, i_m, i_M)
            S = transition(model, m, s, x, M, parms)
            V = [(SVector(M..., el...)-a)./(b.-a) for el in S]
            trembling_hand!(view(Π,tuple(i_m,(Colon() for k in 1:(ndims(Π)-1))...)...), V, w)
        end
    end
    Π0 = (reshape(Π,N,N))

    return Π0 
end

function new_transition(model, dp, x0, exo_grid:: EmptyGrid, endo_grid:: UCGrid; exo=nothing)

    parms = SVector(model.calibration[:parameters]...)

    N_m = 1
    N_s = n_nodes(endo_grid)
    N = N_m*N_s
    Π = zeros(N_s, endo_grid.n...)
    s = nodes(endo_grid)

    a = SVector(endo_grid.min...)
    b = SVector(endo_grid.max...)
    i_m = 1
    x = x0.views[1]
    m = SVector(model.calibration[:exogenous]...)
    if !(exo === nothing)
        m = Dolo.repsvec(exo[1], m)   # z0
    end
    for i_M in 1:n_inodes(dp, i_m)
        M = inode(Point, dp, i_m, i_M)
        if !(exo === nothing)
            M = Dolo.repsvec(exo[2], M)   # z1
        end
        w = iweight(dp, i_m, i_M)
        S = transition(model, m, s, x, M, parms)
        S = [(S[n]-a)./(b-a) for n=1:length(S)]
        trembling_hand!(Π, S, w)
    end

    Π0 = (reshape(Π,N,N))

    return Π0
end

"""
Calculates the new distribution μ à τ = t+1 for a given initial distribution μ at τ = t and a given transition matrix.

# Arguments
* `P::Array{Int64, 2}`: Transition matrix.
* `μ0 ::Vector{Float64}`: Initial distribution μ, of the state (exogenous and endogenous), on the grid at τ = t.

# Returns
* `μ0'*P`: New distribution at τ = t+1 .
"""

function new_distribution(P, μ0)
    return μ0'*P
end
