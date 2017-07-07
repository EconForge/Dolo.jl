@compat abstract type AbstractExogenous end

@compat abstract type AbstractProcess <: AbstractExogenous end
@compat abstract type DiscreteProcess <: AbstractProcess end
@compat abstract type ContinuousProcess <: AbstractProcess end
@compat abstract type IIDExogenous <: ContinuousProcess end

@compat abstract type AbstractDiscretizedProcess <: DiscreteProcess end

###
### Discretized process
###


# date-t grid has a known structure
type DiscretizedProcess{TG<:Grid} <: AbstractDiscretizedProcess
    grid::TG
    integration_nodes::Array{Matrix{Float64},1}
    integration_weights::Array{Vector{Float64},1}
end

n_nodes(dp::DiscretizedProcess) = n_nodes(dp.grid)
node(dp::DiscretizedProcess, i) = node(dp.grid, i)
n_inodes(dp::DiscretizedProcess, i::Int) = size(dp.integration_nodes[i], 1)
inode(dp::DiscretizedProcess, i::Int, j::Int) = dp.integration_nodes[i][j, :]
iweight(dp::DiscretizedProcess, i::Int, j::Int) = dp.integration_weights[i][j]

# date-t grid is unstructured
type DiscreteMarkovProcess <: AbstractDiscretizedProcess
    grid::UnstructuredGrid
    transitions::Matrix{Float64}
    values::Matrix{Float64}
    i0::Int
end

discretize(::Type{DiscreteMarkovProcess}, mp::DiscreteMarkovProcess) = mp
discretize(mp::DiscreteMarkovProcess) = mp

DiscreteMarkovProcess(transitions::Matrix{Float64}, values::Matrix{Float64}) =
    DiscreteMarkovProcess(UnstructuredGrid(values), transitions, values)

DiscreteMarkovProcess(grid::UnstructuredGrid, transitions::Matrix{Float64}, values::Matrix{Float64}) =
    DiscreteMarkovProcess(grid, transitions, values, 1)


n_nodes(dp::DiscreteMarkovProcess) = size(dp.values, 1)
n_inodes(dp::DiscreteMarkovProcess, i::Int) = size(dp.values, 1)
inode(dp::DiscreteMarkovProcess, i::Int, j::Int) = dp.values[j, :]
iweight(dp::DiscreteMarkovProcess, i::Int, j::Int) = dp.transitions[i, j]
node(dp::DiscreteMarkovProcess, i) = dp.values[i, :]
default_index(dp::DiscreteMarkovProcess) = dp.i0


function MarkovProduct(mc1::DiscreteMarkovProcess, mc2::DiscreteMarkovProcess)
    Q = gridmake(mc1.values, mc2.values)
    P = fkron(mc1.transitions, mc2.transitions)
    return DiscreteMarkovProcess(P, Q)
end

function MarkovProduct(mcs::DiscreteMarkovProcess...)
    reduce(MarkovProduct, mcs)
end

# date-t grid is empty

type DiscretizedIIDProcess <: AbstractDiscretizedProcess
    # nodes::Matrix{Float64}
    grid::EmptyGrid
    integration_nodes::Matrix{Float64}
    integration_weights::Vector{Float64}
end

DiscretizedIIDProcess(x, w) = DiscretizedIIDProcess(EmptyGrid(), x, w)

n_nodes(dp::DiscretizedIIDProcess) = 0
n_inodes(dp::DiscretizedIIDProcess, i::Int) = size(dp.integration_nodes, 1)
inode(dp::DiscretizedIIDProcess, i::Int, j::Int) = dp.integration_nodes[j, :]
iweight(dp::DiscretizedIIDProcess, i::Int, j::Int) = dp.integration_weights[j]
node(dip::DiscretizedIIDProcess, i) = zeros(n_inodes(dip, 1))

# Normal law

immutable MvNormal <: IIDExogenous
    mu::Vector{Float64}
    Sigma::Matrix{Float64}
end

MvNormal(Sigma::Matrix{Float64}) = MvNormal(zeros(size(Sigma, 1)), Sigma)
MvNormal(sigma::Float64) = MvNormal(reshape([sigma], 1, 1))

function discretize(mvn::MvNormal)
    n = fill(5, size(mvn.mu))
    x, w = QE.qnwnorm(n, mvn.mu, diagm(vec(mvn.Sigma)))
    DiscretizedIIDProcess(x, w)
end

function Base.rand(mvn::MvNormal, args...)
    dist = Distributions.MvNormal(mvn.mu, mvn.Sigma)
    return rand(dist, args...)
end

function simulate(mvn::MvNormal, N::Integer, T::Int, e0::Vector{Float64}; stochastic::Bool=true)
    dist = QE.MVNSampler(mvn.mu, mvn.Sigma)
    d = length(mvn.mu)
    out = zeros(d, N, T)
    for i in 1:N
        out[:, i, 1] = e0
    end
    if stochastic
        for t in 2:T
            out[:, :, t] = rand(dist, N)
        end
    end
    AxisArray(out, Axis{:V}(1:d), Axis{:N}(1:N), Axis{:T}(1:T))
end

function simulate(mvn::MvNormal, N::Integer, T::Int; stochastic::Bool=true)
    simulate(mvn, N, T, mvn.mu; stochastic=stochastic)
end

function response(mvn::MvNormal, e1::AbstractVector; T::Int=40)
    d = length(mvn.mu)
    out = zeros(d, T)
    out[:, 2] = e1
    return out
end

function simulate(process::DiscreteMarkovProcess, N::Int, T::Int, i0::Int)
    mc_qe = QE.MarkovChain(process.transitions)
    inds = Array{Int}(T, N)
    QE.simulate_indices!(inds, mc_qe, init=i0)
    AxisArray(inds, Axis{:T}(1:T),  Axis{:N}(1:N))
end

function simulate(process::DiscreteMarkovProcess, N::Int, T::Int, m0::AbstractVector{Float64})
    # try to find index in process.values. IF we do, then simulate using
    # the corresponding row index. If we don't, then throw an error
    for i in 1:size(process.values, 1)
        if process.values[i, :] == m0
            return simulate(process, N, T, i)
        end
    end
    error("Couldn't find the vector `m0` in process.values. "*
          "Try passing an integer `i0` instead")
end

function simulate_values(process::DiscreteMarkovProcess, N::Int, T::Int, i0::Int)
    inds = simulate(process, N, T, i0)
    n_values = size(process.values, 2)
    out_values = Array{Float64}(n_values, N, T)
    for n in 1:N, t in 1:T
        out_values[:, n, t] = process.values[inds[t, n], :]
    end
    AxisArray(out_values, Axis{:n}(1:n_values), Axis{:T}(1:T), Axis{:N}(1:N))
end

discretize(dmp::DiscreteMarkovProcess) = dmp


# VAR 1

type VAR1 <: ContinuousProcess
    mu::Array{Float64,1}
    R::Array{Float64,2}
    Sigma::Array{Float64,2}
end

function VAR1(R::Array{Float64,2}, Sigma::Array{Float64,2})
    M = zeros(size(R, 1))
    assert(size(R)==size(Sigma))
    return VAR1(M, R, Sigma)
end

function VAR1(rho::Array{Float64,1}, Sigma::Array{Float64,2})
    R = diagm(rho)
    return VAR1(R, Sigma)
end

function VAR1(rho::Float64, Sigma::Array{Float64,2})
    R = eye(size(Sigma,1))*rho
    return VAR1(R, Sigma)
end

function discretize(var::VAR1)
    d = size(var.R, 1)
    n_states = ones(Int, d)*5
    n_integration = ones(Int, d)*5
    dis = discretize(var, n_states, n_integration)
    return dis
end

function discretize(var::VAR1, n_states::Array{Int,1}, n_integration::Array{Int,1}; n_std::Int=2)
    R = var.R
    M = var.mu
    Sigma = var.Sigma
    S = QE.solve_discrete_lyapunov(R,Sigma)  # asymptotic variance
    sig = diag(S)
    min = var.mu - n_std*sqrt.(sig)
    max = var.mu + n_std*sqrt.(sig)
    grid = CartesianGrid(min,max,n_states)
    # discretize innovations
    x,w = QE.qnwnorm(n_integration, zeros(size(var.Sigma,1)), var.Sigma)
    integration_nodes = [ cat(1,[(M + R*(node(grid, i)-M) + x[j,:])' for j=1:size(x,1)]...) for i in 1:n_nodes(grid)]
    integration_weights = [w for i in 1:n_nodes(grid)]
    return DiscretizedProcess(grid, integration_nodes, integration_weights)
end

discretize(::Type{DiscretizedProcess}, var::VAR1; args...) = discretize(var; args...)

function discretize(::Type{DiscreteMarkovProcess}, var::VAR1; N::Int=3)

    # it would be good to have a special type of VAR1 process
    # which has a scalar autoregressive term
    R =  var.R
    d = size(R, 1)
    ρ = R[1, 1]
    @assert maximum(abs, R-eye(d)*ρ)<1e-16

    sigma = var.Sigma

    if size(var.Sigma, 1) == 1
        mc_qe = QE.rouwenhorst(N, ρ, sqrt(sigma[1]))
        return DiscreteMarkovProcess(mc_qe.p, appenddim(collect(mc_qe.state_values)))
    end


    L = chol(sigma) # sigma = L'*L
    NN = fill(N, d) # default: same number of nodes in each dimension

    components = [QE.rouwenhorst(N, ρ, 1.0) for N in NN]
    mc_components = [DiscreteMarkovProcess(mc.p, appenddim(collect(mc.state_values))) for mc in components]
    mc_prod = MarkovProduct(mc_components...)
    mc_prod.values = mc_prod.values*L'
    return mc_prod
end

"""
```julia
simulate(var::VAR1, N::Int, T::Int, x0::Vector{Float64}
         stochastic::Bool=true, irf::Bool=false, e0::Vector{Float64}=zeros(0))
```

Construct `N` simulated paths of length `T` the vector autoregressive process
`var`, starting each simulated path from `x0`.

If `irf` is `true`, then `e0` must be a vector of the same length as `var.mu`
and `e0` will be used as the initial shocks.

If `stochastic` is `true`, construct a sequence of shocks by drawing `N` paths
of length `T` from the distribution `Normal(0, var.Sigma)` (where `0` is a
vector of zeros of the appropriate length). Otherwise, if `stochastic` is
false, set all shocks to 0 (unless `irf` is true -- see above).

The output is an `AxisArray` contining variables, paths, and time on the 3 axes
"""
function simulate(var::VAR1, N::Int, T::Int, x0::Vector{Float64};
                  stochastic::Bool=true, irf::Bool=false, e0::Vector{Float64}=zeros(0))

    n = size(var.mu, 1)
    XN = zeros(N, T, n)
    E = zeros(N, T, n)

    if stochastic
        dist = QE.MVNSampler(zeros(n), var.Sigma)
        for jj in 1:N
            E[jj, :, :] = rand(dist, T)'
        end
    end

    if irf
        E[:, 1, :] = e0
    end

    # Initial conditions
    for i in 1:N
        XN[i, 1, :] = x0
        for jj in 1:N, ii in 1:T-1
            XN[jj, ii+1, :] = var.mu+var.R*(XN[jj, ii, :]-var.mu)+E[jj, ii, :]
        end
    end
    XN = permutedims(XN, [3, 1, 2])
    AxisArray(XN, Axis{:V}(1:n), Axis{:N}(1:N), Axis{:T}(1:T))
end


function simulate(var::VAR1, N::Int, T::Int; kwargs...)
    return simulate(var, N, T, var.mu; kwargs...)
end

function response(var::VAR1, x0::AbstractVector,
                  e1::AbstractVector; T::Int=40)
    simulate(var, 1, T, x0; stochastic=false, irf=true, e0=e1)[1, :, :]
end

function response(var::VAR1, e1::AbstractVector; T::Int=40)
    simulate(var, 1, T; stochastic=false, irf=true, e0=e1)[1, :, :]
end

function response(var::VAR1, x0::AbstractVector, index_s::Int; T::Int=40)
    e1 = zeros(size(var.mu, 1))
    Impulse = sqrt.(diag(var.Sigma)[index_s])
    e1[index_s] = Impulse
    simulate(var, 1, T, x0; stochastic=false, irf=true, e0=e1)[1, :, :]
end

function response(var::VAR1; T::Int=40)
    e1 = sqrt.(diag(var.Sigma))
    simulate(var, 1, T; stochastic=false, irf=true, e0=e1)[1, :, :]
end

function ErgodDist(var::VAR1, N::Int, T::Int)

    # Simulate for T period, N times
    #  #    - simulate a VAR1 process for T periods N times#
    n = size(var.mu, 1)
    # srand(123) # Setting the seed
    d = MvNormal(var.Sigma)
    VAR_process = simulate(var, N, T)
    E = VAR_process[2]
    XN = VAR_process[1]
    # Computing moments
    Mean_sim = mean(squeeze(mean(VAR_process, 3), 2), 2)
    #Computing the std of simulted Processes (Covariance matrix)
    E1_d = (E[:, :, 1] - repeat(mean(E[:, :, 1], 2), inner=[1, T]))
    Sigma_sim = mean(diag(cov(E1_d[:, :])))  # mean of covariances across simulations

    # Autocorrelation matrix
    X_d = zeros(N, T, 2)
    X_d0 = zeros(N, T-1, 2)
    X_d1 = zeros(N, T-1, 2)
    R_sim = zeros(2, 2)

    for ii in 1:2
        X_d[:, :, ii] = (XN[:, :, 2] - repeat(mean(XN[:, :, ii], 2), inner=[1, T]))
        X_d0[:, :, ii] = X_d[:, 1:end-1, ii]
        X_d1[:, :, ii] = (XN[:, 2:end, ii] - repeat(mean(XN[:, 2:end, ii], 2), inner=[1, T-1]))
        R_sim[ii, ii] = mean(diag( cov(X_d0[:, :, ii], X_d1[:, :, ii] )/sqrt.(var(X_d0[:, :, ii]))/sqrt.(var(X_d1[:, :, ii]))))
        if ii == 2
            R_sim[ii-1, ii] = mean(diag(cov(X_d0[:, :, ii-1], X_d1[:, :, ii])/sqrt.(var(X_d0[:, :, ii-1]))/sqrt.(var(X_d1[:, :, ii]))))
            R_sim[ii, ii-1] = mean(diag(cov(X_d0[:, :, ii], X_d1[:, :, ii-1])/sqrt.(var(X_d0[:, :, ii]))/sqrt.(var(X_d1[:, :, ii-1]))))
        end
    end

    return Mean_sim, Sigma_sim, R_sim
end

#### ProductProcess


type ProductProcess <: AbstractProcess
    process_1::AbstractProcess
    process_2::AbstractProcess
end

function discretize(::Type{DiscreteMarkovProcess}, pp::ProductProcess; opt1=Dict(), opt2=Dict())
    p1 = discretize(DiscreteMarkovProcess, pp.process_1; opt1...)
    p2 = discretize(DiscreteMarkovProcess, pp.process_2; opt2...)
    return MarkovProduct(p1,p2)
end

function discretize(pp::ProductProcess; kwargs...)
    return discretize(DiscreteMarkovProcess, pp; kwargs...)
end


## special processes (could be implemented as types later on)

function DeathProcess(mu::Float64)
    values = [0.0 1.0;]'
    transitions = [(1-mu) mu; 0 1]
    DiscreteMarkovProcess(transitions, values)
end
function PoissonProcess(mu::Float64, K::Int)
    values = (0:K)[:,:]*1.0
    transitions = zeros(K+1, K+1)
    for i=1:K
        transitions[i,i] = 1-mu
        transitions[i,i+1] = mu
    end
    transitions[K+1,K+1] = 1
    DiscreteMarkovProcess(transitions, values)
end
function AgingProcess(mu::Float64, K::Int)
    values = zeros(K+1,2)
    values[:,1] = 0:K
    values[1,2] = 1
    transitions = zeros(K+1, K+1)
    transitions[1,1] = 1
    for i=2:K
        transitions[i,i+1] = (1-mu)
        transitions[i,1] = mu
    end
    transitions[end,1] = 1
    dp = DiscreteMarkovProcess(transitions, values)
    dp.i0 = 2
    dp
end


# compatibility names
@compat const AR1 = VAR1
@compat const MarkovChain = DiscreteMarkovProcess
@compat const Normal = MvNormal
@compat const GDP = DiscretizedProcess
