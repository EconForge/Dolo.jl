abstract type AbstractExogenous end

abstract type EmptyProcess <: AbstractExogenous end
abstract type AbstractProcess <: AbstractExogenous end
abstract type DiscreteProcess <: AbstractProcess end
abstract type ContinuousProcess <: AbstractProcess end
abstract type IIDExogenous <: ContinuousProcess end

abstract type AbstractDiscretizedProcess <: DiscreteProcess end

# this is a bit crude but not performance critical, for now
function node(::Type{Point}, dp::DiscreteProcess, i::Int)
    v = node(dp, i)
    SVector{length(v)}(v)
end
function inode(::Type{Point}, dp::DiscreteProcess, i::Int, j::Int)
    v = inode(dp, i, j)
    SVector{length(v)}(v)
end

struct ConstantProcess <: ContinuousProcess
    μ::Array{Float64,1}
end

ConstantProcess(;μ=zeros(1)) = ConstantProcess(μ)
ConstantProcess(;μ=zeros(1)) = ConstantProcess(μ)



###
### Discretized process
###

# date-t grid has a known structure
mutable struct GDP{TG} <: AbstractDiscretizedProcess
    grid::TG
    integration_nodes::Array{Matrix{Float64},1}
    integration_weights::Array{Vector{Float64},1}
end

n_nodes(dp::GDP) = n_nodes(dp.grid)
node(dp::GDP, i) = node(dp.grid, i)
n_inodes(dp::GDP, i::Int) = size(dp.integration_nodes[i], 1)
inode(dp::GDP, i::Int, j::Int) = dp.integration_nodes[i][j, :]
iweight(dp::GDP, i::Int, j::Int) = dp.integration_weights[i][j]


# function node(::Type{Point}, dp::DiscreteProcess, i::Int) where n_x
#     SVector(dp.grid.nodes[i]...)
# end
# function inode(::Type{Point}, dp::DiscreteProcess, i::Int, j::Int) where n_x
#     v = inode(dp, i, j)
#     SVector(dp.grid.integration_nodes[i][j,:]...)
# end

function node(::Type{Point{n_x}}, dp::GDP, i::Int) where n_x
    SVector{n_x}(dp.grid.nodes[i]...)
end
function inode(::Type{Point{n_x}}, dp::GDP, i::Int, j::Int) where n_x
    SVector{n_x, Float64}(dp.integration_nodes[i][j,:]...)
end


function Product(gdp1::GDP, gdp2::GDP)
  In = [ gridmake(gdp1.integration_nodes[i], gdp2.integration_nodes[j] ) for i = 1:n_nodes(gdp1)  for j = 1:n_nodes(gdp2) ]
  Iw =[kron(gdp1.integration_weights[i], gdp2.integration_weights[j] ) for i = 1:n_nodes(gdp1)  for j = 1:n_nodes(gdp2) ]
  
  i_grid  = ProductGrid(gdp1.grid, gdp2.grid)

  return GDP(i_grid, In, Iw)
end


# date-t grid is unstructured
mutable struct DiscreteMarkovProcess <: AbstractDiscretizedProcess
    grid::UnstructuredGrid
    transitions::Matrix{Float64}
    values::Matrix{Float64}
    i0::Int
end

ndims(dmp::DiscreteMarkovProcess) = ndims(dmp.grid)

discretize(::Type{DiscreteMarkovProcess}, mp::DiscreteMarkovProcess) = mp
discretize(mp::DiscreteMarkovProcess) = mp

DiscreteMarkovProcess(transitions::Matrix{Float64}, values::Matrix{Float64}) =
    DiscreteMarkovProcess(UnstructuredGrid{size(values,2)}(values), transitions, values)

DiscreteMarkovProcess(grid::UnstructuredGrid, transitions::Matrix{Float64}, values::Matrix{Float64}) =
    DiscreteMarkovProcess(grid, transitions, values, 1)


n_nodes(dp::DiscreteMarkovProcess) = size(dp.values, 1)
n_inodes(dp::DiscreteMarkovProcess, i::Int) = size(dp.values, 1)
iweight(dp::DiscreteMarkovProcess, i::Int, j::Int) = dp.transitions[i, j]
default_index(dp::DiscreteMarkovProcess) = dp.i0
inode(dp::DiscreteMarkovProcess, i::Int, j::Int) = dp.values[j, :]
node(dp::DiscreteMarkovProcess, i) = dp.values[i, :]
inode(::Type{Point{d}}, dp::DiscreteMarkovProcess, i::Int, j::Int) where d = SVector{d}(dp.values[j, :]...)
node(::Type{Point{d}}, dp::DiscreteMarkovProcess, i) where d = SVector{d}( dp.values[i, :] ...)




function MarkovProduct(mc1::DiscreteMarkovProcess, mc2::DiscreteMarkovProcess)
    Q = gridmake(mc1.values, mc2.values)
    P = fkron(mc1.transitions, mc2.transitions)
    return DiscreteMarkovProcess(P, Q)
end


function MarkovProduct(mcs::DiscreteMarkovProcess...)
    reduce(MarkovProduct, mcs)
end

# date-t grid is empty

mutable struct DiscretizedIIDProcess <: AbstractDiscretizedProcess
    # nodes::Matrix{Float64}
    grid::EmptyGrid
    integration_nodes::Matrix{Float64}
    integration_weights::Vector{Float64}
end

DiscretizedIIDProcess(x, w) = DiscretizedIIDProcess(EmptyGrid{size(x,2)}(), x, w)

n_nodes(dp::DiscretizedIIDProcess) = 0
n_inodes(dp::DiscretizedIIDProcess, i::Int) = size(dp.integration_nodes, 1)
inode(dp::DiscretizedIIDProcess, i::Int, j::Int) = dp.integration_nodes[j, :]
iweight(dp::DiscretizedIIDProcess, i::Int, j::Int) = dp.integration_weights[j]
node(dip::DiscretizedIIDProcess, i::Int) = fill(NaN, n_inodes(dip, 1))

node(::Type{Point{d}}, dip::DiscretizedIIDProcess, i::Int) where d = fill(0, SVector{d,Float64})
# node(::Type{Point{d}}, dip::DiscretizedIIDProcess, i::Int) where d = fill(NaN, SVector{d,Float64}) # TODO: this should be the correct version !!

inode(::Type{Point{d}}, dip::DiscretizedIIDProcess, i::Int, j::Int) where d = SVector{d}(dip.integration_nodes[j, :]...)

# Normal law

struct MvNormal <: IIDExogenous
    μ::Vector{Float64}
    Σ::Matrix{Float64}
end

MvNormal(Σ::Matrix{Float64}) = MvNormal(zeros(size(Σ, 1)), Σ)
MvNormal(σ::Float64) = MvNormal(reshape([σ^2], 1, 1))

Normal(;Σ=zeros(1,1)) = MvNormal(Σ)

UNormal(;σ=0.0) = MvNormal(reshape([σ^2], 1, 1))


function discretize(mvn::MvNormal; n=5::Union{Int, Vector{Int}})
    x, w = QE.qnwnorm(n, mvn.μ, mvn.Σ)
    DiscretizedIIDProcess(x, w)
end


function Base.rand(mvn::MvNormal, args...)
    dist = Distributions.MvNormal(mvn.μ, mvn.Σ)
    return rand(dist, args...)
end

function simulate(mvn::MvNormal, e0::Vector{Float64}; N::Integer, T::Int, stochastic::Bool=true)
    dist = QE.MVNSampler(mvn.μ, mvn.Σ)
    d = length(mvn.μ)
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

function simulate(mvn::MvNormal, N::Integer, T::Int, e0::Vector{Float64}; stochastic::Bool=true)
    simulate(mvn, e0; N=N, T=T, stochastic=stochastic)
end

function simulate(mvn::MvNormal; N::Integer, T::Int, stochastic::Bool=true)
    simulate(mvn, mvn.μ; N=N, T=T, stochastic=stochastic)
end

function simulate(mvn::MvNormal, N::Integer, T::Int; stochastic::Bool=true)
    simulate(mvn, N, T, mvn.μ; stochastic=stochastic)
end

function response(mvn::MvNormal, e1::AbstractVector; T::Int=40)
    d = length(mvn.μ)
    out = zeros(d, T)
    out[:, 2] = e1
    return out
end

function simulate(process::DiscreteMarkovProcess; N::Int, T::Int, i0::Int)
    mc_qe = QE.MarkovChain(process.transitions)
    inds = Array{Int}(undef, T, N)
    QE.simulate_indices!(inds, mc_qe, init=i0)
    AxisArray(inds, Axis{:T}(1:T),  Axis{:N}(1:N))
end

function simulate(process::DiscreteMarkovProcess, N::Int, T::Int, i0::Int)
    simulate(process; N=N, T=T, i0=i0)
end

function simulate(process::DiscreteMarkovProcess, m0::AbstractVector{Float64}; N::Int, T::Int)
    # try to find index in process.values. IF we do, then simulate using
    # the corresponding row index. If we don't, then throw an error
    for i in 1:size(process.values, 1)
        if process.values[i, :] == m0
            return simulate(process; N=N, T=T, i0=i)
        end
    end
    error("Couldn't find the vector `m0` in process.values. "*
          "Try passing an integer `i0` instead")
end

function simulate(process::DiscreteMarkovProcess, N::Int, T::Int, m0::AbstractVector{Float64})
    simulate(process, m0; N=N, T=T)
end

function simulate_values(process::DiscreteMarkovProcess; N::Int, T::Int, i0::Int)
    inds = simulate(process; N=N, T=T, i0=i0)
    n_values = size(process.values, 2)
    out_values = Array{Float64}(n_values, N, T)
    for n in 1:N, t in 1:T
        out_values[:, n, t] = process.values[inds[t, n], :]
    end
    AxisArray(out_values, Axis{:n}(1:n_values), Axis{:T}(1:T), Axis{:N}(1:N))
end

function simulate_values(process::DiscreteMarkovProcess, N::Int, T::Int, i0::Int)
    simulate_values(process; N=N, T=T, i0=i0)
end

# VAR 1

mutable struct VAR1 <: ContinuousProcess
    μ::Array{Float64,1}
    R::Array{Float64,2}
    Σ::Array{Float64,2}
end

VAR1(;ρ::Float64=0.0, Σ=ones(1,1)) = VAR1(ρ, Σ)

function VAR1(R::Array{Float64,2}, Σ::Array{Float64,2})
    M = zeros(size(R, 1))
    @assert size(R)==size(Σ)
    return VAR1(M, R, Σ)
end

function VAR1(ρ::Array{Float64,1}, Σ::Array{Float64,2})
    R = diagm(ρ)
    return VAR1(R, Σ)
end

function VAR1(ρ::Float64, Σ::Array{Float64,2})
    p = size(Σ,1)
    R = Matrix(ρ*I, p, p)
    return VAR1(R, Σ)
end

function discretize(::Type{GDP}, var::VAR1; n::Union{Int, Vector{Int}}=5, n_i::Union{Int, Vector{Int}}=5,  n_std::Int=2)
    R = var.R
    M = var.μ
    Σ = var.Σ
    S = QE.solve_discrete_lyapunov(R,Σ)  # asymptotic variance
    sig = diag(S)
    min = var.μ - n_std*sqrt.(sig)
    max = var.μ + n_std*sqrt.(sig)
    if n isa Int
        N = fill(n, length(min))  ::Vector{Int}
    else
        N = n ::Vector{Int}
    end
    grid = UCGrid{length(min)}(min,max,N)
    # discretize innovations
    x,w = QE.qnwnorm(n_i, zeros(size(var.Σ,1)), var.Σ)
    integration_nodes = [ cat([(M + R*(node(grid, i)-M) + x[j,:])' for j=1:size(x,1)]...; dims=1) for i in 1:n_nodes(grid)]
    integration_weights = [w for i in 1:n_nodes(grid)]
    return GDP(grid, integration_nodes, integration_weights)
end

discretize(var::VAR1; args...) = discretize(DiscreteMarkovProcess, var; args...)


function discretize(::Type{DiscreteMarkovProcess}, var::VAR1; n::Union{Int, Vector{Int}}=3)

    # it would be good to have a special type of VAR1 process
    # which has a scalar autoregressive term
    R =  var.R
    d = size(R, 1)
    ρ = R[1, 1]
    @assert maximum(abs, R-Matrix(ρ*I,d,d))<1e-16

    σ = var.Σ

    if size(var.Σ, 1) == 1
        mc_qe = QE.rouwenhorst(n, ρ, sqrt(σ[1]))
        return DiscreteMarkovProcess(mc_qe.p, appenddim(collect(mc_qe.state_values)))
    end


    L = chol(σ) # σ = L'*L
    if n isa Int
        NN = fill(n, d) # default: same number of nodes in each dimension
    else
        NN = n
    end

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

If `irf` is `true`, then `e0` must be a vector of the same length as `var.μ`
and `e0` will be used as the initial shocks.

If `stochastic` is `true`, construct a sequence of shocks by drawing `N` paths
of length `T` from the distribution `Normal(0, var.Σ)` (where `0` is a
vector of zeros of the appropriate length). Otherwise, if `stochastic` is
false, set all shocks to 0 (unless `irf` is true -- see above).

The output is an `AxisArray` contining variables, paths, and time on the 3 axes
"""
function simulate(var::VAR1, x0::Vector{Float64}; N::Int, T::Int,
    stochastic::Bool=true, irf::Bool=false, e0::Vector{Float64}=zeros(0))

    n = size(var.μ, 1)
    XN = zeros(n, N, T)
    E = zeros(n, N, T)

    if stochastic
        dist = QE.MVNSampler(zeros(n), var.Σ)
        for jj in 1:N
            E[:, jj, :] = rand(dist, T)
        end
    end

    if irf
        E[:, :, 1] = repeat(e0,N,1)
    end

    # Initial conditions
    for i in 1:N
        XN[:, i, 1] = x0
        for  ii in 1:T-1
            XN[:, i, ii+1] = var.μ+var.R*(XN[:, i, ii]-var.μ)+E[:, i, ii]
        end
    end
    AxisArray(XN, Axis{:V}(1:n), Axis{:N}(1:N), Axis{:T}(1:T))
end

function simulate(var::VAR1, N::Int, T::Int, x0::Vector{Float64};
    stochastic::Bool=true, irf::Bool=false, e0::Vector{Float64}=zeros(0))

    simulate(var, x0; N=N, T=T, stochastic=stochastic, irf=irf, e0=e0)
end

function simulate(var::VAR1; N::Int, T::Int, kwargs...)
    return simulate(var, var.μ; N=N, T=T, kwargs...)
end

function simulate(var::VAR1, N::Int, T::Int; kwargs...)
    return simulate(var; N=N, T=T, kwargs...)
end

function response(var::VAR1, x0::AbstractVector,
                  e1::AbstractVector; T::Int=40)
    simulate(var, x0; N=1, T=T, stochastic=false, irf=true, e0=e1)[1, :, :]
end

function response(var::VAR1, e1::AbstractVector; T::Int=40)
    simulate(var; N=1, T=T, stochastic=false, irf=true, e0=e1)[1, :, :]
end

function response(var::VAR1, x0::AbstractVector, index_s::Int; T::Int=40)
    e1 = zeros(size(var.μ, 1))
    Impulse = sqrt.(diag(var.Σ)[index_s])
    e1[index_s] = Impulse
    simulate(var, x0; N=1, T=T, stochastic=false, irf=true, e0=e1)[1, :, :]
end

function response(var::VAR1; T::Int=40)
    e1 = sqrt.(diag(var.Σ))
    simulate(var; N=1, T=T, stochastic=false, irf=true, e0=e1)[1, :, :]
end

function ErgodDist(var::VAR1, N::Int, T::Int)

    # Simulate for T period, N times
    #  #    - simulate a VAR1 process for T periods N times#
    n = size(var.μ, 1)
    # srand(123) # Setting the seed
    d = MvNormal(var.Σ)
    VAR_process = simulate(var, N, T)
    E = VAR_process[2]
    XN = VAR_process[1]
    # Computing moments
    Mean_sim = mean(squeeze(mean(VAR_process, 3), 2), 2)
    #Computing the std of simulted Processes (Covariance matrix)
    E1_d = (E[:, :, 1] - repeat(mean(E[:, :, 1], 2), inner=[1, T]))
    Σ_sim = mean(diag(cov(E1_d[:, :])))  # mean of covariances across simulations

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

    return Mean_sim, Σ_sim, R_sim
end

#### ProductProcess


mutable struct ProductProcess{P1<:AbstractProcess,P2<:AbstractProcess} <: AbstractProcess
    process_1::P1
    process_2::P2
end

ProductProcess(p) = p

function discretize(pp::ProductProcess{ConstantProcess, <:IIDExogenous}; opt=Dict())
    diidp = discretize(pp.process_2)
    inodes = diidp.integration_nodes
    iit = hcat(
        vcat([pp.process_1.μ' for i=1:size(inodes,1)]...), 
        inodes
    )
    return DiscretizedIIDProcess(diidp.grid, iit, diidp.integration_weights)
end

function discretize(::Type{DiscreteMarkovProcess}, pp::ProductProcess; opt1=Dict(), opt2=Dict())
    p1 = discretize(DiscreteMarkovProcess, pp.process_1; opt1...)
    p2 = discretize(DiscreteMarkovProcess, pp.process_2; opt2...)
    return MarkovProduct(p1,p2)
end


function discretize(::Type{GDP}, pp::ProductProcess; opt1=Dict(), opt2=Dict())
  p1 = discretize(GDP, pp.process_1; opt1...)
  p2 = discretize(GDP, pp.process_2; opt2...)
  return Product(p1,p2)
end


function discretize(pp::ProductProcess; kwargs...)
    return discretize(DiscreteMarkovProcess, pp; kwargs...)
end


## special processes (could be implemented as types later on)

function DeathProcess(μ::Float64)
    values = [0.0 1.0;]'
    transitions = [(1-μ) μ; 0 1]
    DiscreteMarkovProcess(transitions, values)
end
function PoissonProcess(μ::Float64, K::Int)
    values = (0:K)[:,:]*1.0
    transitions = zeros(K+1, K+1)
    for i=1:K
        transitions[i,i] = 1-μ
        transitions[i,i+1] = μ
    end
    transitions[K+1,K+1] = 1
    DiscreteMarkovProcess(transitions, values)
end
function AgingProcess(μ::Float64, K::Int)
    values = zeros(K+1,2)
    values[:,1] = 0:K
    values[1,2] = 1
    transitions = zeros(K+1, K+1)
    transitions[1,1] = 1
    for i=2:K
        transitions[i,i+1] = (1-μ)
        transitions[i,1] = μ
    end
    transitions[end,1] = 1
    dp = DiscreteMarkovProcess(transitions, values)
    dp.i0 = 2
    dp
end


get_integration_nodes(dprocess::Dolo.AbstractDiscretizedProcess, i::Int)= [(iweight(dprocess,i,j), inode(dprocess,i,j), j) for j in 1:n_inodes(dprocess,i) if iweight(dprocess,i,j)!=0]

get_integration_nodes(::typeof(Point), dprocess::Dolo.AbstractDiscretizedProcess, i::Int)=[(iweight(dprocess,i,j), inode(Point,dprocess,i,j), j) for j in 1:n_inodes(dprocess,i) if iweight(dprocess,i,j)!=0]


# type unstable
#
# get_integration_nodes(dprocess::Dolo.AbstractDiscretizedProcess, i::Int)=Iterators.filter( x -> (x[1]!=0), ((iweight(dprocess,i,j), inode(dprocess,i,j), j) for j in 1:n_inodes(dprocess,i)) )
#
# get_integration_nodes(::typeof(Point), dprocess::Dolo.AbstractDiscretizedProcess, i::Int)=Iterators.filter( x -> (x[1]!=0), ((iweight(dprocess,i,j), inode(Point,dprocess,i,j), j) for j in 1:n_inodes(dprocess,i)) )


# compatibility names
const AR1 = VAR1
const MarkovChain = DiscreteMarkovProcess
# const DiscretizedProcess

MarkovChain(;transitions=ones(1,1), values=[range(1,size(transitions,1))...]) = MarkovChain(transitions, values)





#### Domains

get_domain(mv::MvNormal) = EmptyDomain{length(mv.μ)}()

function get_domain(dmp::DiscreteMarkovProcess)
    points = dmp.grid.nodes
    d = ndims(dmp)
    return DiscreteDomain{d}(points)
end

function get_domain(cp::ConstantProcess)
    d = length(cp.μ) 
    min  = fill(-Inf, d) ### not absolutely clear this is the right default
    max  = fill(Inf, d)
    return CartesianDomain(min, max)
end

function get_domain(iid::DiscretizedIIDProcess)
    return EmptyDomain()
end

function get_domain(var::VAR1)
    d = length(var.μ)
    return CartesianDomain(fill(-Inf, d), fill(Inf, d))
end

function get_domain(pp::ProductProcess)
    dom1 = get_domain(pp.process_1)
    dom2 = get_domain(pp.process_2)
    return ProductDomain(dom1, dom2)
end


# function discretize(pp::ProductProcess{ConstantProcess, VAR1}; opts...)
#     cp = pp.process_1
#     gdp = discretize(pp.process_2; opts...)
#     v = SVector(cp.μ...)
#     grid = [SVector(v..., e...) for e in gdp.grid.nodes]
#     integration_nodes = [cat(fill(cp.μ, size(e,1)), e; dims=2) for e in gdp.integration_nodes]
#     return DiscretizedProcess(grid, integration_nodes, gdp.integration_weights)

# end

function discretize(::Type{DiscreteMarkovProcess},pp::ProductProcess{ConstantProcess, VAR1}; opts...)
    mc = discretize(DiscreteMarkovProcess, pp.process_2; opts...)
    cp = pp.process_1
    d1 = length(cp.μ)
    d2 = length(node(mc.grid,1))
    v = SVector(cp.μ...)
    grid = [SVector(v..., e...) for e in mc.grid.nodes]
    transitions = mc.transitions
    values = copy(from_LOP(grid))
    return DiscreteMarkovProcess(transitions, values)
end





function discretize(::Type{GDP}, p::ConstantProcess)
    μ = p.μ
    grid = PointGrid(μ)
    integration_nodes = [repeat(μ',1)]
    integration_weights = [[1.0]]
    return GDP(grid, integration_nodes, integration_weights)
end

