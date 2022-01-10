function (G::distG)(μ0::AbstractVector{Float64}, x0::MSM{Point{n_x}}; exo =nothing, diff=false) where n_x


    if !diff
        P = transition_matrix(G.model, G.dprocess, x0, G.grid; exo=exo, diff=false)
        μ1 = P'*μ0
        return μ1
    end

    P, P_x = transition_matrix(G.model, G.dprocess, x0, G.grid; exo=exo, diff=true)

    μ1 = P'μ0
    
    M = length(μ0)
    N = length(x0.data)*length(x0.data[1])

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

    ∂G_∂μ = LinearMap( μ -> P'*μ, M, M)
    ∂G_∂x = LinearMap(fun_x, M, N)
    return μ1, ∂G_∂μ, ∂G_∂x

end

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
            dm_dz = Dolo.repsvec((@SVector ones(n_z1)),m*0)
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
                dM_dz = Dolo.repsvec((@SVector ones(n_2)),M*0)
            end
            w = iweight(dp, i_m, i_M)
            if diff
                S, S_z1, S_x, S_z2 = transition(model, Val{(0,1,3,4)}, m, s, x, M, parms)
            else
                S = transition(model, m, s, x, M, parms)
            end
            S = [(S[n]-a)./(b-a) for n=1:length(S)]
            ind_s = tuple((Colon() for k in 1:(ndims(Π)-3))...)
            Π_view = view(Π,:,i_m,ind_s..., i_MM)
            if !diff
                trembling_hand!(Π_view, S, w)
            else
                S_x = [( 1.0 ./(b-a)) .* S_x[n] for n=1:length(S)]
                S_z1 = [( 1.0 ./(b-a)) .* S_z1[n] for n=1:length(S)] * dm_dz
                S_z2 = [( 1.0 ./(b-a)) .* S_z2[n] for n=1:length(S)] * dM_dz

                dΠ_view_x = view(dΠ_x,:,i_m,ind_s...,i_MM)
                dΠ_view_z1 = view(dΠ_z1,:,i_m,ind_s...,i_MM)
                dΠ_view_z2 = view(dΠ_z2,:,i_m,ind_s...,i_MM)

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


function transition_matrix(model, sol; diff=false)
    x0 = Dolo.MSM([sol.dr(i, sol.dr.grid_endo.nodes) for i=1:max(1,Dolo.n_nodes(sol.dr.grid_exo))])
    grid = ProductGrid(sol.dr.grid_exo, sol.dr.grid_endo)
    Dolo.transition_matrix(model, sol.dprocess, x0, grid; diff=diff);
end

function transition_matrix(G::distG; dp=G.dprocess, x0=G.x0, grid=G.grid, exo=nothing, diff=false)

    return transition_matrix(G.model, dp, x0, grid; exo=exo, diff=diff)
end


function trembling_foot!(Π, dΠ_x, dΠ_z1, dΠ_z2, S::Vector{Point{d}}, S_x::Vector{SMatrix{d,n_x,Float64,_}}, S_z1, S_z2, w::Float64) where d where n_x where _
    
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
        
        @assert d==1
        
        λ_vec =  (SVector( -1. /δ[1], 1. / δ[1]), )
        A = outer(λ_vec...)

        rhs_dΠ_x = outer2(A, Sn_x[1,:])
        rhs_dΠ_x = [ -1. /δ[1].*Sn_x[1,:] , 1. / δ[1].*Sn_x[1,:]]
        dΠ_x[indexes_to_be_modified...] .+= w*rhs_dΠ_x

        rhs_dΠ_z1 = outer2(A, Sn_z1[1,:])
        rhs_dΠ_z1 = [ -1. /δ[1].*Sn_z1[1,:] , 1. / δ[1].*Sn_z1[1,:]]
        dΠ_z1[indexes_to_be_modified...] .+= w*rhs_dΠ_z1

        rhs_dΠ_z2 = outer2(A, Sn_z2[1,:])
        rhs_dΠ_z2 = [ -1. /δ[1].*Sn_z2[1,:] , 1. / δ[1].*Sn_z2[1,:]]
        dΠ_z2[indexes_to_be_modified...] .+= w*rhs_dΠ_z2
        
    end

end
