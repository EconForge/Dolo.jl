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



function (G::distG)(μ0::AbstractVector{Float64}, x0::MSM{Point{n_x}}; exo =nothing, diff=false) where n_x

    if !diff
        P = transition_matrix(G.model, G.dprocess, x0, G.grid; exo=exo, diff=false)
        μ1 = P'*μ0
        return μ1
    end

    P, P_x, P_z1, P_z2 = transition_matrix(G.model, G.dprocess, x0, G.grid; exo=exo, diff=true)

    #dev
    #P_z2 = @SMatrix [reshape([P_z2[:,15*i+j] for j=1:15, i=0:14]',1,225)[i][j] for j=1:225, i=1:225]

    μ1 = P'μ0
    
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

    print(" typeof Pz2 ",typeof(P_z2)," ", "size of Pz2", size(P_z2)," ")
    print(" typeof P ",typeof(P)," ", "size of P", size(P)," ")

    function fun_z2(dz2)
        P_dz2 = [(P_z2[i,j]'*dz2)  for i=1:size(P,1), j=1:size(P,2)]
        print("typeof Pdz2 ",typeof(P_dz2)," size of P_dz2 ", size(P_dz2)," ")
        d_μ = P_dz2'*μ0
        return d_μ
    end


    ∂G_∂μ = LinearMap( μ -> P'*μ, M, M)
    ∂G_∂x = LinearMap(fun_x, M, N)
    ∂G_∂z1 = LinearMap(fun_z1, M, N_z1)
    ∂G_∂z2 = LinearMap(fun_z2, M, N_z2)

    return μ1, ∂G_∂μ, ∂G_∂x, ∂G_∂z1, ∂G_∂z2

end

function (G::distG)(μ0::AbstractVector{Float64}, x0::AbstractVector{Float64}; exo=nothing, diff=false)
    
    n_x = length(G.x0.data[1])
    x = MSM(copy(reinterpret(SVector{n_x, Float64}, x0)), G.x0.sizes)
    return G(μ0, x; exo=exo, diff=diff)

end


"""
Updates A.
# Arguments
* `A`: the transition matrix that will be updated.
* `x` : vector of controls.
* `w::Float64` : vector of weights.
* `exo_grid` : exogenous grid that will be used to determine the type of rescaling to do.
* `a` : SVector containing the minimum values on the UCGrids
* `b` : SVector containing the maximum values on the UCGrids
# Optional 
* `M` : future node considered when the exogenous grid is a UCGrid to rescale x
# Modifies
* `A` : the updated transition matrix 
"""
function trembling_hand_rescaled!(A, x, w::Float64, exo_grid, a, b; M=0)
    if typeof(exo_grid) == Dolo.UnstructuredGrid{ndims(exo_grid)}
        x = [(x[n]-a)./(b-a) for n=1:length(x)]
        trembling_hand!(A, x, w)
    elseif typeof(exo_grid) == Dolo.UCGrid{ndims(exo_grid)}
        V = [(SVector(M..., el...)-a)./(b.-a) for el in x]
        trembling_hand!(A, V, w)
    else

        x = [(x[n]-a)./(b-a) for n=1:length(x)]
        trembling_hand!(A, x, w)
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

function transition_matrix(model, dp, x0::MSM{<:SVector{n_x}}, grid; exo=nothing, diff=false) where n_x

    # the restriction here is that dp is a descrete_process

    parms = SVector(model.calibration[:parameters]...)
    
    exo_grid = grid.exo
    endo_grid = grid.endo

    N_m = max(1, n_nodes(exo_grid))
    N_s = n_nodes(endo_grid)
    N = N_m*N_s
    Π = zeros(N_s, N_m, endo_grid.n..., N_m)
    if diff
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
                #println("M : ",M," ")
            end
            w = iweight(dp, i_m, i_M)
            if diff
                S, S_z1, S_x, S_z2 = transition(model, Val{(0,1,3,4)}, m, s, x, M, parms)
            else
                S = transition(model, m, s, x, M, parms)
            end
            print("S_z2 = ", S_z2," ")
            S = [(S[n]-a)./(b-a) for n=1:length(S)]
            ind_s = tuple((Colon() for k in 1:(ndims(Π)-3))...)
            Π_view = view(Π,:,i_m,ind_s..., i_MM)
            if !diff
                trembling_hand!(Π_view, S, w)
            else
                #print("type of Sz1 0 : ",typeof(S_z1)," length of Sz1 0 : ", length(S_z1), " ") # Vector{SMatrix{2, 1, Float64, 2}} of length 225
                #print("type of Sx 0 : ",typeof(S_x)," length of Sx 0 : ", length(S_x), " ") # Vector{SMatrix{2, 2, Float64, 4}} of length 225

                S_x = [( 1.0 ./(b-a)) .* S_x[n] for n=1:length(S)]
                S_z1 = [( 1.0 ./(b-a)) .* S_z1[n] for n=1:length(S)] 
                S_z2 = [( 1.0 ./(b-a)) .* S_z2[n] for n=1:length(S)] 

                #print("type of Sz1 1 : ",typeof(S_z1)," length of Sz1 1 : ", length(S_z1), " ") # Vector{SMatrix{2, 1, Float64, 2}} of length 225
                #print("type of Sx 1 : ",typeof(S_x)," length of Sx 1 : ", length(S_x), " ") # Vector{SMatrix{2, 2, Float64, 4}} of length 225

                S_z1 = [M[:,1:length(exo[1])] for M in S_z1]
                S_z2 = [M[:,1:length(exo[2])] for M in S_z2] 

                #print("type of Sz1 2 : ",typeof(S_z1)," length of Sz1 2 : ", length(S_z1), " size of Sz1[1] 2",size(S_z1[1])," ") #  Vector{Matrix{Float64}} of length 225 containing matrices of size (2,1) for rbc_iid

                dΠ_view_x = view(dΠ_x,:,i_m,ind_s...,i_MM)
                dΠ_view_z1 = view(dΠ_z1,:,i_m,ind_s...,i_MM)
                dΠ_view_z2 = view(dΠ_z2,:,i_m,ind_s...,i_MM)

                print("dΠdΠ_view_z2_z2 : ", dΠ_view_z2," ")

                trembling_foot!(Π_view, dΠ_view_x, dΠ_view_z1, dΠ_view_z2, S, S_x, S_z1, S_z2, w)
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


function transition_matrix(model, sol; diff=false, exo=nothing)
    x0 = Dolo.MSM([sol.dr(i, sol.dr.grid_endo.nodes) for i=1:max(1,Dolo.n_nodes(sol.dr.grid_exo))])
    grid = ProductGrid(sol.dr.grid_exo, sol.dr.grid_endo)
    Dolo.transition_matrix(model, sol.dprocess, x0, grid; diff=diff, exo=nothing);
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


function trembling_foot!(Π, dΠ_x, dΠ_z1, dΠ_z2, S::Vector{Point{d}}, S_x::Vector{SMatrix{d,n_x,Float64,_}}, S_z1, S_z2, w::Float64) where d where n_x where _
    
    @assert ndims(Π) == d+1
    shape_Π = size(Π)
    grid_dimension = d
    δ =  SVector{d,Float64}(1.0./(shape_Π[1+i]-1) for i in 1:d )
    N = shape_Π[1]
    
    #print("typeof Snx : ", typeof(S_x[1]), " ") # SMatrix{2, 2, Float64, 4} 
    #print("typeof Snz2 : ", typeof(S_z2[1]), "size of Snz2 : ",size(S_z2[1])," ") # Matrix{Float64} (2,1)

    for n in 1:N

        Sn = S[n]
        Sn_x = S_x[n]
        Sn_z1 = S_z1[n]
        Sn_z2 = S_z2[n]

        inbound = (0. .<= Sn) .& ( Sn .<= 1.0)

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
    
        #print("d = ",d," ") # 2
        for k=1:d
            
            λ_vec =  tuple( (i==k ? SVector( -1. /δ[k] * inbound[k], 1. / δ[k] * inbound[k]) : (SVector((1-λn[i]),λn[i])) for i in 1:d)... )
            A = outer(λ_vec...)

            rhs_dΠ_x = outer2(A, Sn_x[k,:])
            dΠ_x[indexes_to_be_modified...] .+= w*rhs_dΠ_x


            #print("size of rhs_dΠ_x : ", size(rhs_dΠ_x)," typeof rhs_dΠ_x", typeof(rhs_dΠ_x)," ") #SizedMatrix{2, 2, SVector{2, Float64}, 2, Matrix{SVector{2, Float64}}}

            rhs_dΠ_z1 = outer2(A, Sn_z1[k,:])
            dΠ_z1[indexes_to_be_modified...] .+= w*rhs_dΠ_z1

            rhs_dΠ_z2 = outer2(A, Sn_z2[k,:])
            dΠ_z2[indexes_to_be_modified...] .+= w*rhs_dΠ_z2
            #print("size of rhs_dΠ_z2 : ", size(rhs_dΠ_z2)," typeof rhs_dΠ_z2", typeof(rhs_dΠ_z2)," ") #SizedMatrix{2, 2, Vector{Float64}, 2, Matrix{Vector{Float64}}} 

        end
        # dΠ_z2_view = view(dΠ_z2,n,1,:,:,1) 
        # dΠ_z2_view = view(dΠ_z2,n,1,:,:,1)'
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

