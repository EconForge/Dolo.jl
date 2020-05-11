import Dolo: Model, time_iteration, discretize
using LinearAlgebra, Plots


model = Model("examples/models/rbc_mc.yaml")


sol = time_iteration(model)


import Dolo: UnstructuredGrid, CartesianGrid, transition, n_nodes, n_inodes, inode, iweight, nodes, node, Point

dp = discretize(model.exogenous)


function ergodic_distribution(model, dr, exo_grid:: UnstructuredGrid, endo_grid:: CartesianGrid, dp)
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
            trembling_hand!(view(Π,i_m, :, i_M, :), S, w)
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


function trembling_hand!(A, x, w)
    N,n0 = size(A)
    δ0 = 1.0./(n0-1.0)
    for n in 1:N
        print("hey")
        x0 = x[n][1]
        x0 = min.(max.(x0, 0.0),1.0)
        q0 = div.(x0, δ0)
        q0 = max.(0, q0)
        q0 = min.(q0, n0-2)
        λ0 = (x0./δ0-q0) # ∈[0,1[ by construction
        q0_ = round.(Int,q0) + 1
        println("n:", n, ":", q0_, ":", λ0, ":", w)
        A[n, q0_]   += (1-λ0)*w
        A[n, q0_+1] += λ0*w
    end
end


Π, μ = ergodic_distribution(model, dr, dr.grid_exo, dr.grid_endo, dp)

sum(reshape(Π, 40,40), dims=2)

μ = reshape(μ, 2, 20)

plot(μ[1,:])


#
# function trembling_hand(A, x, w)
#
#     N,n0,n1 = A.shape
#     δ0 = 1/(n0-1)
#     δ1 = 1/(n1-1)
#
#     for n in 1:N
#
#         x0 = x[n,0]
#         x0 = min(max(x0, 0),1)
#         q0 = np.floor_divide(x0, δ0)
#         q0 = max(0, q0)
#         q0 = min(q0, n0-2)
#
#         x1 = x[n,1]
#         x1 = min(max(x1, 0),1)
#         q1 = np.floor_divide(x1, δ1)
#         q1 = max(0, q1)
#         q1 = min(q1, n1-2)
#
#         λ0 = (x0-q0*δ0)/δ0 # ∈[0,1[ by construction
#         q0_ = int(q0)
#
#         λ1 = (x1-q1*δ1)/δ1 # ∈[0,1[ by construction
#         q1_ = int(q1)
#
#         A[n, q0_, q1_]   += (1-λ0)*(1-λ1)*w
#         A[n, q0_, q1_+1]   += (1-λ0)*(λ1)*w
#         A[n, q0_+1, q1_]   += (λ0)*(1-λ1)*w
#         A[n, q0_+1, q1_+1]   += (λ0)*(λ1)*w
#
#     end
# end
#




dr(1,s[1])

ss0

s[1]
