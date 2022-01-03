# function label_density(μ, symbols,  grid)
#     exo_names = 
#     endo_names = symbols[:states]
#     sc = scales(grid_endo)
#     d = Dict(endo_names[i]=>sc[i] for i=1:length(grid_endo.n))
#     return AxisArray(μ; d...)
# end


function ergodic_distribution(model, sol) where T where Q
    G = distG(model, sol)
    μ = ergodic_distribution(G)
    return μ
    # return label_density(μ, model.symbols,  sol.dr.grid_exo, sol.dr.grid_endo)
end


# represents the way that distributions are discretized
struct distG{ID, n_m, n_s, n_x, Gx, Ge}

    model::AbstractModel

    grid::ProductGrid{Gx, Ge}
    dprocess

    s0::ListOfPoints{n_s}
    x0::MSM{SVector{n_x, Float64}}
    μ0::AbstractVector{Float64}

end

function distG(model, sol)
    ID = id(model)
    grid_endo = sol.dr.grid_endo
    grid_exo = sol.dr.grid_exo
    grid = ProductGrid(grid_exo, grid_endo)
    dprocess = sol.dprocess
    s0 = grid_endo.nodes
    x0 = MSM([sol.dr(i, s0) for i=1:max(1,n_nodes(grid_exo))])
    N = length(grid)
    μ0 = ones(N)/N
    n_m = ndims(grid_exo)
    n_s = ndims(grid_endo)
    n_x = length(x0.data[1])
    Gx = typeof(grid_exo)
    Ge = typeof(grid_endo)
    return distG{ID, n_m, n_s, n_x, Gx, Ge}(model, grid, dprocess, s0, x0, μ0)
end


function ergodic_distribution(G::distG; x0=G.x0, kwargs...)
    P = transition_matrix(G.model, G.dprocess, x0, G.grid; kwargs...)
    M = [P'-I; ones(1, size(P,1))]
    R = zeros(size(P,1)+1)
    R[end] = 1
    μ = M \ R  
    return μ
end



function (G::distG)(μ0::AbstractVector{Float64}, x0::MSM{Point{n_x}}; diff=false) where n_x


    if !diff
        P = transition_matrix(G.model, G.dprocess, x0, G.grid; exo=nothing, diff=false)
        μ1 = P'*μ0
        return μ1
    end

    P, P_x = transition_matrix(G.model, G.dprocess, x0, G.grid; exo=nothing, diff=true)

    μ1 = P'*μ0
    
    M = length(μ0)
    N = length(x0.data)*length(x0.data[1])

    function fun_x(dx::AbstractVector{Float64})
        d_x = MSM(copy(reinterpret(SVector{n_x, Float64}, dx)), x0.sizes)
        P_dx = [(P_x[i,j]'*d_x.data[i])  for i=1:size(P,1), j=1:size(P,2)]
        return P_dx'*μ0
    end

    ∂G_∂μ = LinearMap( μ -> P'*μ, M, M)
    ∂G_∂x = LinearMap(fun_x, M, N)
    return μ1, ∂G_∂μ, ∂G_∂x

end

function (G::distG)(μ0::AbstractVector{Float64}, x0::AbstractVector{Float64}; diff=false)
    
    n_x = length(G.x0.data[1])
    x = MSM(copy(reinterpret(SVector{n_x, Float64}, x0)), G.x0.sizes)
    return G(μ0, x; diff=diff)

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

function transition_matrix(model, dp, x0::MSM{<:SVector{n_x}}, grid; exo=nothing, diff=false) where n_x

    # the restriction here is that dp is a descrete_process

    parms = SVector(model.calibration[:parameters]...)
    
    exo_grid = grid.exo
    endo_grid = grid.endo

    N_m = max(1, n_nodes(exo_grid))
    N_s = n_nodes(endo_grid)
    N = N_m*N_s
    Π = zeros(N_m, N_s, N_m, endo_grid.n...)
    if diff
        dΠ = zeros(SVector{n_x, Float64}, N_m, N_s, N_m, endo_grid.n...)
    end
    s = nodes(endo_grid)
    a = SVector(endo_grid.min...)
    b = SVector(endo_grid.max...)
    for i_m in 1:max(1, n_nodes(exo_grid))
        x = x0.views[i_m]
        m = node(exo_grid, i_m)
        if !(exo === nothing)
            m = Dolo.repsvec(exo[1], m)   # z0
        end
        for i_M in 1:n_inodes(dp, i_m)
            
            if isa(grid.exo, EmptyGrid)
                i_MM = 1 # used to index transition matrix
            else
                i_MM = i_M
            end
            M = inode(Point, dp, i_m, i_M)
            if !(exo === nothing)
                M = Dolo.repsvec(exo[2], M)   # z1
            end
            w = iweight(dp, i_m, i_M)
            println(i_m, " ; " ,i_M, " ; " , i_MM, " ; " , w)
            if diff
                S, S_x = transition(model, Val{(0,3)}, m, s, x, M, parms)
            else
                S = transition(model, m, s, x, M, parms)
            end
            S = [(S[n]-a)./(b-a) for n=1:length(S)]
            Π_view = view(Π,tuple(i_m,:,i_MM,(Colon() for k in 1:(ndims(Π)-3))...)...)
            if !diff
                trembling_hand!(Π_view, S, w)
            else
                S_x = [( 1.0 ./(b-a)) .* S_x[n] for n=1:length(S)]
                dΠ_view = view(dΠ,tuple(i_m,:,i_MM,(Colon() for k in 1:(ndims(dΠ)-3))...)...)
                trembling_foot!(Π_view, dΠ_view, S, S_x, w)
            end
        end
    end

    Π0 = (reshape(Π,N,N))
    if !diff
        return Π0
    else
        dΠ0 = reshape(dΠ,N,N)
        return Π0, dΠ0
    end

end


function transition_matrix(model, sol; diff=false)
    x0 = Dolo.MSM([sol.dr(i, sol.dr.grid_endo.nodes) for i=1:max(1,Dolo.n_nodes(sol.dr.grid_exo))])
    grid = ProductGrid(sol.dr.grid_exo, sol.dr.grid_endo)
    Dolo.transition_matrix(model, sol.dprocess, x0, grid; diff=diff);
end

function transition_matrix(G::distG; dp=G.dprocess, x0=G.x0, grid=G.grid, exo=nothing, diff=false)

    return transition_matrix(G.model, dp, x0, grid; exo=exo, diff=diff)
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
        qn_ = round.(Int,qn) .+ 1
        
        λn_weight_vector = tuple( (SVector((1-λn[i]),λn[i]) for i in 1:d)... )

        indexes_to_be_modified = tuple(n, UnitRange.(qn_,qn_.+1)...)

        # Filling transition matrix
        rhs = outer(λn_weight_vector...)
        A[indexes_to_be_modified...] .+= w*rhs
        
    end

end


function trembling_foot!(Π, dΠ, S::Vector{Point{d}}, S_x::Vector{SMatrix{d,n_x,Float64,_}}, w::Float64) where d where n_x where _
    
    @assert ndims(Π) == d+1
    shape_Π = size(Π)
    grid_dimension = d
    δ =  SVector{d,Float64}(1.0./(shape_Π[1+i]-1) for i in 1:d )
    N = shape_Π[1]

    for n in 1:N

        Sn = S[n]
        Sn_x = S_x[n]

        Sn = min.(max.(Sn, 0.0),1.0)
        qn = div.(Sn, δ)
        qn = max.(0, qn)
        qn = min.(qn, shape_Π[2:d+1].-2)
        λn = (Sn./δ.-qn) # ∈[0,1[ by construction
        qn_ = round.(Int,qn) .+ 1
        
        λn_weight_vector_Π = tuple( (SVector((1-λn[i]),λn[i]) for i in 1:d)... )

        # # # Filling transition matrix
        rhs_Π = outer(λn_weight_vector_Π...)

        indexes_to_be_modified = tuple(n, UnitRange.(qn_,qn_.+1)...)

        Π[indexes_to_be_modified...] .+= w.*rhs_Π

        for k=1:d
            λ_vec =  tuple( (i==k ? SVector( -1. /(size(Π, k+1)-1), 1. /(size(Π, k+1)-1) ) : (SVector((1-λn[i]),λn[i])) for i in 1:d)... )
            A = outer(λ_vec...)
            rhs_dΠ = outer2(A, Sn_x[k,:])
            dΠ[indexes_to_be_modified...] .+= w*rhs_dΠ
        end
        
    end

end



"""
Computes the outer product.

TODO: define the outer product

# Argument
* `λn_weight_vector::Vararg{Point{2}}`: tuple of Point{2} to be multiplied by outer product

# Returns
* the outer product
"""
function outer(λn_weight_vector::Vararg{SVector{2}})
    return [prod(e) for e in Iterators.product(λn_weight_vector...)]
end

# TODO: define and document
outer2(A, x) = [A[i]*x for i in CartesianIndices(A)]

