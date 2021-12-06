using StaticArrays
using Dolo
using Dolo: Point
pwd()
filename = "./examples/models/rbc.yaml"
model = Model(filename)
sol = time_iteration(model)
dr = sol.dr
dr([1.],[2.1])
grid, dp = Dolo.discretize(model)
Dolo.discretize(Dolo.DiscreteMarkovProcess, model.exogenous)
F = Dolo.Euler(model)
F.x0
Dolo.discretize(model.exogenous)
grid
typeof(dr.grid_endo)

transition(model, Dolo.node(grid.exo, 1), Dolo.nodes(grid.endo), dr(1, Dolo.nodes(grid.endo)),Dolo.inode(Point, dp, 1, 1),SVector(model.calibration[:parameters]...))

function new_distribution_v1(sol, μ0, x0, dr, parms)
    endo_grid, exo_grid = dr.grid_endo, dr.grid_exo
    dp = sol.dprocess
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

function new_distribution_v2(sol, μ0, x0, exo_grid:: UnstructuredGrid, endo_grid:: UCGrid, parms)

    dp = sol.dprocess

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

function new_distribution_v2(sol, μ0, x0, exo_grid:: UCGrid, endo_grid:: UCGrid, parms)

    dp = sol.dprocess

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
            trembling_hand!(view(Π,tuple(i_m,:,i_M,(Colon() for k in 1:(ndims(Π)-1))...)...), V, w)
        end
    end
    Π0 = (reshape(Π,N,N))

    return Π0 * μ0
end

function new_distribution_v2(sol, μ0, x0, exo_grid:: EmptyGrid, endo_grid:: UCGrid, parms)

    dp = sol.dprocess

    N_m = 1
    N_s = n_nodes(endo_grid)
    N = N_m*N_s
    Π = zeros(N_s, endo_grid.n...)
    s = nodes(endo_grid)

    a = SVector(endo_grid.min...)
    b = SVector(endo_grid.max...)
    i_m = 1
    x = dr(s)
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

C = zeros(2, 3, 2, 3)
xc = [SVector{3,Float64}(3.4,4.7,1.1),SVector{3,Float64}(2.6,4.2,1.0)]
wc =1.0
my_trembling_hand_v3!(C,xc,wc)
print(C)
print(3)

B = Array{Float64,3}(undef, 2, 2, 2)
xb = [SVector{2,Float64}(3.4,4.7),SVector{2,Float64}(2.6,4.2)]
wb =1.0
my_trembling_hand_v3!(B,xb,wb)
B
view(B,tuple(1,(Colon() for k in 1:2)...)...)
