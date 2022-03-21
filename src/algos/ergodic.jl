# function label_density(μ, symbols,  grid)
#     exo_names = 
#     endo_names = symbols[:states]
#     sc = scales(grid_endo)
#     d = Dict(endo_names[i]=>sc[i] for i=1:length(grid_endo.n))
#     return AxisArray(μ; d...)
# end

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


"""
Computes the ergodic distribution of the states (endo & exo).
# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `sol`: solution of the model obtained after a time iteration algorithm
# Optional Argument
* `smooth::boolean`: Indicates whether, when we discretize the  transition results on the grid, we want to smooth the initial linear ponderation of nodes or not.
# Returns
* `μ1::distG`: the ergodic distribution
"""
function ergodic_distribution(model, sol; smooth=true) where T where Q
    G = distG(model, sol)
    μ = ergodic_distribution(G; smooth=smooth)
    return μ
    # return label_density(μ, model.symbols,  sol.dr.grid_exo, sol.dr.grid_endo)
end

function ergodic_distribution(G::distG; x0=G.x0, kwargs...)
    P = transition_matrix(G.model, G.dprocess, x0, G.grid; kwargs...)
    M = [P'-I; ones(1, size(P,1))]
    R = zeros(size(P,1)+1)
    R[end] = 1
    μ = M \ R  
    return μ
end


"""
Uses the current distribution of the states (endogeneous & exogeneous) and the controls to compute the next distribution of the states on the grid. Exogeneous parameters can be included
as well as the possibility to differentiate this function with respect to the states distribution, the vector of controls and the exogeneous parameters.
# Arguments
* `μ0::AbstractVector{Float64}`: current distribution of the states (endo & exo)
* `x0::MSM{Point{n_x}}`: vector of controls
# Optional Arguments
* `exo`: nothing or (z0, z1)
* `diff::boolean`: Indicates whether we want to evaluate the derivative of G wrt μ, x and possibly z1 and z2
* `smooth::boolean`: Indicates whether, when we discretize the  transition results on the grid, we want to smooth the initial linear ponderation of nodes or not.
# Returns
* `μ1::distG`: It is the new distribution of the states (endogeneous & exogeneous)
# Optionally returns also
* `∂G_∂μ::LinearMaps.FunctionMap{Float64}`: It is the jacobian of G wrt μ.
* `∂G_∂x::LinearMaps.FunctionMap{Float64}`: It is the jacobian of G wrt x.
* `∂G_∂z1::LinearMaps.FunctionMap{Float64}`: It is the jacobian of G wrt z1.
* `∂G_∂z2::LinearMaps.FunctionMap{Float64}`: It is the jacobian of G wrt z2.
"""
function (G::distG)(μ0::AbstractVector{Float64}, x0::MSM{Point{n_x}}; exo =nothing, diff=false, smooth=true) where n_x

    if !diff
        P = transition_matrix(G.model, G.dprocess, x0, G.grid; exo=exo, diff=false)
        μ1 = P'*μ0
        return μ1
    end

    P, P_x, P_z1, P_z2 = transition_matrix(G.model, G.dprocess, x0, G.grid; exo=exo, diff=true, smooth=smooth)


    μ1 = P'μ0 #new distribution of the states
    
    M = length(μ0)
    N = length(x0.data)*length(x0.data[1])
    N_z1 = length(exo[1])
    N_z2 = length(exo[2])

    function fun_x(dx::AbstractVector{Float64})
        d_x = MSM(copy(reinterpret(SVector{n_x, Float64}, dx)), x0.sizes)
        N = size(P,2)
        d_μ = [ sum([μ0[k]*(P_x[k,i]'*d_x.data[k]) for k=1:N]) for i=1:N ]
        # alternate calculation
        P_dx = [(P_x[i,j]'*d_x.data[i])  for i=1:size(P,1), j=1:size(P,2)]
        d_μμ = P_dx'*μ0
        @assert (maximum(abs, d_μ-d_μμ)<1e-10)
        return d_μ
    end

    function fun_z1(dz1) 
        P_dz1 = [(P_z1[i,j]'*dz1)  for i=1:size(P,1), j=1:size(P,2)]
        d_μ = P_dz1'*μ0
        return d_μ
    end

    function fun_z2(dz2)
        P_dz2 = [(P_z2[i,j]'*dz2)  for i=1:size(P,1), j=1:size(P,2)]
        d_μ = P_dz2'*μ0
        return d_μ
    end


    ∂G_∂μ = LinearMap( μ -> P'*μ, M, M)
    ∂G_∂x = LinearMap(fun_x, M, N)
    ∂G_∂z1 = LinearMap(fun_z1, M, N_z1)
    ∂G_∂z2 = LinearMap(fun_z2, M, N_z2)

    return μ1, ∂G_∂μ, ∂G_∂x, ∂G_∂z1, ∂G_∂z2

end

function (G::distG)(μ0::AbstractVector{Float64}, x0::AbstractVector{Float64}; exo=nothing, diff=false, smooth=true)
    
    n_x = length(G.x0.data[1])
    x = MSM(copy(reinterpret(SVector{n_x, Float64}, x0)), G.x0.sizes)
    return G(μ0, x; exo=exo, diff=diff, smooth=smooth)

end



"""
Calculates the new transition matrix for a given model, a given discretized exogenous process, given control values (x0) and given grids (exogenous and endogenous).
# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `dprocess::`: Discretized exogenous process.
* `x0::Dolo.MSM{SVector{2, Float64}}`: Initial control values.
* `grid`: grid
# Optional arguments
* `exo`: nothing or (z0, z1)
* `diff`: true or false. Indicates whether we want to evaluate the derivative of G wrt μ, x and possibly z1 and z2
* `smooth`: true or false. Indicates whether, when we discretize the  transition results on the grid, we want to smooth the initial linear ponderation of nodes or not.
# Returns
* `Π0::`: New transition matrix.
"""

function transition_matrix(model, dp, x0::MSM{<:SVector{n_x}}, grid; exo=nothing, diff=false, smooth = true) where n_x

    # the restriction here is that dp is a descrete_process

    parms = SVector(model.calibration[:parameters]...)
    
    exo_grid = grid.exo
    endo_grid = grid.endo

    N_m = max(1, n_nodes(exo_grid))
    N_s = n_nodes(endo_grid) 
    N = N_m*N_s
    Π = zeros(N_s, N_m, endo_grid.n..., N_m) #initialization of the transition matrix
    if diff #initialization of its derivatives
        dΠ_x = zeros(SVector{n_x, Float64}, N_s, N_m, endo_grid.n..., N_m)
        if !(exo === nothing)
            n_z1 = length(exo[1])
            n_z2 = length(exo[2])
            dΠ_z1 = zeros(SVector{n_z1, Float64}, N_s, N_m, endo_grid.n..., N_m)
            dΠ_z2 = zeros(SVector{n_z2, Float64}, N_s, N_m, endo_grid.n..., N_m)
        end
    end
    s = nodes(endo_grid)
    a = SVector(endo_grid.min...)
    b = SVector(endo_grid.max...)
    for i_m in 1:max(1, n_nodes(exo_grid))
        x = x0.views[i_m]
        m = node(exo_grid, i_m)
        if !(exo === nothing)
            m = Dolo.repsvec(exo[1], m)   # completes z0 by elements of m until reaching a length equal to max(length(m),length(z0))
        end
        for i_M in 1:n_inodes(dp, i_m)
            
            if isa(grid.exo, EmptyGrid)
                i_MM = 1 # used to index transition matrix
            else
                i_MM = i_M
            end
            M = inode(Point, dp, i_m, i_M) # future node
            if !(exo === nothing)
                M = Dolo.repsvec(exo[2], M)   # z1
            end
            w = iweight(dp, i_m, i_M)
            if diff
                S, S_z1, S_x, S_z2 = transition(model, Val{(0,1,3,4)}, m, s, x, M, parms) # computes the next vector of states and its derivatives wrt z1, x and z2. 
                # The points obtained are located between points of the endogenous grid and a discretization needs therefore to be done to adjust the transition matrix.
            else
                S = transition(model, m, s, x, M, parms)
            end
            S = [(S[n]-a)./(b-a) for n=1:length(S)]
            ind_s = tuple((Colon() for k in 1:(ndims(Π)-3))...)
            Π_view = view(Π,:,i_m,ind_s..., i_MM)
            if !diff
                trembling_hand!(Π_view, S, w; smooth=smooth)
            else
                
                S_x = [( 1.0 ./(b-a)) .* S_x[n] for n=1:length(S)]
                S_z1 = [( 1.0 ./(b-a)) .* S_z1[n] for n=1:length(S)] 
                S_z2 = [( 1.0 ./(b-a)) .* S_z2[n] for n=1:length(S)] 

                S_z1 = [M[:,1:length(exo[1])] for M in S_z1]
                S_z2 = [M[:,1:length(exo[2])] for M in S_z2] 

               
                dΠ_view_x = view(dΠ_x,:,i_m,ind_s...,i_MM)
                dΠ_view_z1 = view(dΠ_z1,:,i_m,ind_s...,i_MM)
                dΠ_view_z2 = view(dΠ_z2,:,i_m,ind_s...,i_MM)

                trembling_foot!(Π_view, dΠ_view_x, dΠ_view_z1, dΠ_view_z2, S, S_x, S_z1, S_z2, w; smooth = smooth)
            end
        end
    end

    Π0 = (reshape(Π,N,N))
    if !diff
        return Π0
    else
        dΠ_0_x = reshape(dΠ_x,N,N)
        dΠ_0_z1 = reshape(dΠ_z1,N,N)
        dΠ_0_z2 = reshape(dΠ_z2,N,N)
        return Π0, dΠ_0_x, dΠ_0_z1, dΠ_0_z2
    end

end


function transition_matrix(model, sol; diff=false, exo=nothing, smooth=true)
    x0 = Dolo.MSM([sol.dr(i, sol.dr.grid_endo.nodes) for i=1:max(1,Dolo.n_nodes(sol.dr.grid_exo))])
    grid = ProductGrid(sol.dr.grid_exo, sol.dr.grid_endo)
    Dolo.transition_matrix(model, sol.dprocess, x0, grid; diff=diff, exo=nothing, smooth=smooth);
end

function transition_matrix(G::distG; dp=G.dprocess, x0=G.x0, grid=G.grid, exo=nothing, diff=false, smooth=true)

    return transition_matrix(G.model, dp, x0, grid; exo=exo, diff=diff,smooth=smooth)
end

"""
Evaluate a polynomial that allows to smooth the repartition done by trembling_foot.
# Arguments
* `x`: a float number.
# Returns
* the value of 2x^3 - 3x^2 +1 
"""
function smoothing(x::SVector{2})
    return 2. .* x .^3 .- 3. .* x.^2 .+1.
end

"""
Evaluate the derivative of a polynomial that allows to smooth the repartition done by trembling_foot.
# Arguments
* `x`: a float number.
# Returns
* the value of the derivative of x |—> 2x^3 - 3x^2 +1 
"""
function ∂smoothing(x::SVector{2})
    return 6. .* x .^2 .- 6. .* x 
end

"""
Updates A.
# Arguments
* `A`: the transition matrix that will be updated.
* `x::Vector{Point{d}}` : vector of controls.
* `w::Float64` : vector of weights.
# Optional Argument
* `smooth::boolean` : indicates whether the discretization is made with a smooth ponderation or a linear one
# Modifies
* `A` : the updated transition matrix 
"""
function trembling_hand!(A, x::Vector{Point{d}}, w::Float64; smooth=true) where d
    
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
        
        if smooth==false
            λn_weight_vector =  tuple( (SVector((1-λn[i]),λn[i]) for i in 1:d)... )
        else
            λn_weight_vector = tuple( (smoothing(SVector((λn[i]),1-λn[i])) for i in 1:d)... )
        end
        
        indexes_to_be_modified = tuple(n, UnitRange.(qn_,qn_.+1)...)

        # Filling transition matrix
        rhs = outer(λn_weight_vector...)
        A[indexes_to_be_modified...] .+= w*rhs
        
    end

end


"""
Updates Π, the transition matrix, and its derivatives with respect to x, z1 and z2. This is done using the previous Π, dΠ_x, dΠ_z1 and dΠ_z2 and the derivatives
of the transition function wrt x, z1 and z2.
# Arguments
* `Π`: the transition matrix that will be updated.
* `dΠ_x` : its derivative wrt x
* `dΠ_z1` : its derivative wrt z1
* `dΠ_z2` : its derivative wrt z2
* `S::Vector{Point{d}}` : Vector of SVectors of states (endo & exo) with d the dimension of the state space and the size of the vector corresponding to the size of the grid
* `S_x` : Derivative of S wrt x
* `S_z1` : Derivative of S wrt z1
* `S_z2` : Derivative of S wrt z2
* `w::Float64` : weight
# Optional Argument
* `smooth::boolean` : indicates whether the discretization is made with a smooth ponderation or a linear one
# Modifies
* `Π`: the transition matrix that will be updated.
* `dΠ_x`: its derivative wrt x
* `dΠ_z1`: its derivative wrt z1
* `dΠ_z2`: its derivative wrt z2
"""

function trembling_foot!(Π, dΠ_x, dΠ_z1, dΠ_z2, S::Vector{Point{d}}, S_x::Vector{SMatrix{d,n_x,Float64,_}}, S_z1, S_z2, w::Float64; smooth=true) where d where n_x where _
    
    @assert ndims(Π) == d+1
    shape_Π = size(Π)
    grid_dimension = d
    δ =  SVector{d,Float64}(1.0./(shape_Π[1+i]-1) for i in 1:d )
    N = shape_Π[1]
    

    for n in 1:N

        Sn = S[n]
        Sn_x = S_x[n]
        Sn_z1 = S_z1[n]
        Sn_z2 = S_z2[n]

        inbound = (0. .<= Sn) .& ( Sn .<= 1.0) # detects the points that have fallen outside of the grid

        Sn = min.(max.(Sn, 0.0),1.0)
        qn = div.(Sn, δ)
        qn = max.(0, qn)
        qn = min.(qn, shape_Π[2:d+1].-2)
        λn = (Sn./δ.-qn) # ∈[0,1[ by construction
        qn_ = round.(Int,qn) .+ 1
        
        if smooth==false
            λn_weight_vector_Π = tuple( (SVector((1-λn[i]),λn[i]) for i in 1:d)... ) # linear ponderation
        else
            λn_weight_vector_Π = tuple( (smoothing(SVector((λn[i]),1-λn[i])) for i in 1:d)... ) # smoothed ponderation —> allows to have continuous derivatives 
        end

        # # # Filling transition matrix
        rhs_Π = outer(λn_weight_vector_Π...)

        indexes_to_be_modified = tuple(n, UnitRange.(qn_,qn_.+1)...)

        Π[indexes_to_be_modified...] .+= w.*rhs_Π
        
        # # # Computing the associated derivatives
        for k=1:d
            
            if smooth==false
                λ_vec = tuple( (i==k ? SVector( -1. /δ[k] * inbound[k], 1. / δ[k] * inbound[k]) : (SVector((1-λn[i]),λn[i])) for i in 1:d)... ) 
            else
                λ_vec = tuple( (i==k ? ∂smoothing(SVector((λn[i]),1-λn[i])) .* SVector{2,Float64}(1.,-1.) ./δ[k] .* inbound[k] : (smoothing(SVector((λn[i]),1-λn[i]))) for i in 1:d)... ) 
            end
            
            A = outer(λ_vec...)

            rhs_dΠ_x = outer2(A, Sn_x[k,:])
            dΠ_x[indexes_to_be_modified...] .+= w*rhs_dΠ_x

            rhs_dΠ_z1 = outer2(A, Sn_z1[k,:])
            dΠ_z1[indexes_to_be_modified...] .+= w*rhs_dΠ_z1

            rhs_dΠ_z2 = outer2(A, Sn_z2[k,:])
            dΠ_z2[indexes_to_be_modified...] .+= w*rhs_dΠ_z2
           
        end
    end

end



"""
Computes the outer product.
# Argument
* `λn_weight_vector::Vararg{Point{2}}`: tuple of Point{2} to be multiplied by outer product
# Returns
* the outer product
"""
function outer(λn_weight_vector::Vararg{SVector{2}})
    return [prod(e) for e in Iterators.product(λn_weight_vector...)]
end


"""
Computes an outer product between a vector and a matrix and returns a vector of matrices.
# Arguments
* `A` : SVector of any size
* `x`: matrix of any size
# Returns
* a sized vector obtained by making the outer product of A and x and dividing the matrix obtained in a vector of matrices for which the value of the A[i] term changes
"""
outer2(A, x) = [A[i]*x for i in CartesianIndices(A)]

