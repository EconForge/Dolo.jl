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



function new_distribution(model, sol, μ0, x0, exo_grid:: UnstructuredGrid, endo_grid:: UCGrid)

    dp = sol.dprocess
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
        for i_M in 1:n_inodes(dp, i_m)
            M = inode(Point, dp, i_m, i_M)
            w = iweight(dp, i_m, i_M)
            S = transition(model, m, s, x, M, parms)
            S = [(S[n]-a)./(b-a) for n=1:length(S)]
            trembling_hand!(view(Π,tuple(i_m,:,i_M,(Colon() for k in 1:(ndims(Π)-3))...)...), S, w)
        end
    end
    Π0 = (reshape(Π,N,N))

    return Π0 * μ0
end

function new_distribution(model, sol, μ0, x0, exo_grid:: UCGrid, endo_grid:: UCGrid)

    dp = sol.dprocess
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
        for i_M in 1:n_inodes(dp, i_m)
            M = inode(Point, dp, i_m, i_M)
            w = iweight(dp, i_m, i_M)
            S = transition(model, m, s, x, M, parms)
            V = [(SVector(M..., el...)-a)./(b.-a) for el in S]
            trembling_hand!(view(Π,tuple(i_m,(Colon() for k in 1:(ndims(Π)-1))...)...), V, w)
        end
    end
    Π0 = (reshape(Π,N,N))

    return Π0 * μ0
end

function new_distribution(model, sol, μ0, x0, exo_grid:: EmptyGrid, endo_grid:: UCGrid)

    dp = sol.dprocess
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
    for i_M in 1:n_inodes(dp, i_m)
        M = inode(Point, dp, i_m, i_M)
        w = iweight(dp, i_m, i_M)
        S = transition(model, m, s, x, M, parms)
        S = [(S[n]-a)./(b-a) for n=1:length(S)]
        trembling_hand!(Π, S, w)
    end

    Π0 = (reshape(Π,N,N))

    return Π0 * μ0
end




function trembling_hand!(A::AbstractArray{Float64,2}, x, w)
    N,n0 = size(A)
    δ0 = 1.0./(n0-1.0)
    for n in 1:N
        x0 = x[n][1]
        x0 = min.(max.(x0, 0.0),1.0)
        q0 = div.(x0, δ0)
        q0 = max.(0, q0)
        q0 = min.(q0, n0-2)
        λ0 = (x0./δ0-q0) # ∈[0,1[ by construction
        q0_ = round.(Int,q0) + 1
        A[n, q0_]   += (1-λ0)*w
        A[n, q0_+1] += λ0*w
    end
end

function trembling_hand!(A::AbstractArray{Float64,3}, x, w)
    N,n0,n1 = size(A)
    δ0 = 1.0./(n0-1.0)
    δ1 = 1.0./(n1-1.0)
    for n in 1:N
        x0 = x[n][1]
        x0 = min.(max.(x0, 0.0),1.0)
        q0 = div.(x0, δ0)
        q0 = max.(0, q0)
        q0 = min.(q0, n0-2)
        λ0 = (x0./δ0-q0) # ∈[0,1[ by construction
        q0_ = round.(Int,q0) + 1

        x1 = x[n][2]
        x1 = min.(max.(x1, 0.0),1.0)
        q1 = div.(x1, δ1)
        q1 = max.(0, q1)
        q1 = min.(q1, n1-2)
        λ1 = (x1./δ1-q1) # ∈[0,1[ by construction
        q1_ = round.(Int,q1) + 1

        A[n, q0_ ,  q1_] += (1-λ0)*(1-λ1)*w
        A[n, q0_+1, q1_] += λ0*(1-λ1)*w
        A[n, q0_, q1_+1] += (1-λ0)*λ1*w
        A[n, q0_+1, q1_+1] += λ0*λ1*w
    end
end
