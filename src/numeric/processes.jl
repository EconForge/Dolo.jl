import QuantEcon
const QE = QuantEcon
import Distributions:rand
import Distributions

abstract AbstractExogenous

abstract AbstractProcess <: AbstractExogenous
abstract DiscreteProcess <: AbstractProcess
abstract ContinuousProcess <: AbstractProcess
abstract IIDExogenous <: AbstractProcess

abstract AbstractDiscretizedProcess # <: DiscreteProcess

###
### Discretized process
###

# date-t grid has a known structure
type DiscretizedProcess <: AbstractDiscretizedProcess
    grid::Grid
    integration_nodes::Array{Matrix{Float64},1}
    integration_weights::Array{Vector{Float64},1}
end

n_nodes(dp::DiscretizedProcess) = n_nodes(dp.grid)
node(dp::DiscretizedProcess, i) = node(dp.grid,i)
n_inodes(dp::DiscretizedProcess, i::Int) = size(dp.integration_nodes[i], 1)
inode(dp::DiscretizedProcess, i::Int, j::Int) = dp.integration_nodes[i][j, :]
iweight(dp::DiscretizedProcess, i::Int, j::Int) = dp.integration_weights[i][j]

# date-t grid is unstructured
type DiscreteMarkovProcess <: AbstractDiscretizedProcess
    transitions::Matrix{Float64}
    values::Matrix{Float64}
end

n_nodes(dp::DiscreteMarkovProcess) = size(dp.values, 1)
n_inodes(dp::DiscreteMarkovProcess, i::Int) = size(dp.values, 1)
inode(dp::DiscreteMarkovProcess, i::Int, j::Int) = dp.values[j, :]
iweight(dp::DiscreteMarkovProcess, i::Int, j::Int) = dp.transitions[i,j]
node(dp::DiscreteMarkovProcess, i) = dp.values[i, :]

# function discretize(dmp::DiscreteMarkovProcess)
#     nodes = dmp.values
#     n_nodes = size(nodes, 1)
#     integration_nodes = [nodes for i=1:n_nodes]
#     integration_weights = [dmp.transitions[i, :] for i=1:n_nodes]
#     return DiscretizedProcess(nodes, integration_nodes, integration_weights)
# end
# type DiscreteMarkovProcess <: AbstractDiscretizedProcess
#     nodes::Matrix{Float64}
#     integration_nodes::Array{Matrix{Float64}}
#     integration_weights::Array{Vector{Float64}}
# end
#
# n_nodes(dp::DiscretizedProcess) = size(dp.nodes, 1)
# n_inodes(dp::DiscretizedProcess, i::Int) = size(dp.integration_nodes[i], 1)
# inodes(dp::DiscretizedProcess, i::Int, j::Int) = dp.integration_nodes[i][j, :]
# iweights(dp::DiscretizedProcess, i::Int, j::Int) = dp.integration_weights[i][j]
# node(dp::DiscretizedProcess, i) = dp.nodes[i, :]

# date-t grid is empty
type DiscretizedIIDProcess <: AbstractDiscretizedProcess
    # nodes::Matrix{Float64}
    integration_nodes::Matrix{Float64}
    integration_weights::Vector{Float64}
end

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
    x, w = QE.qnwnorm(n, mvn.mu, mvn.Sigma)
    DiscretizedIIDProcess(x, w)
end

function rand(mvn::MvNormal, args...)
    dist = Distributions.MvNormal(mvn.mu, mvn.Sigma)
    return rand(dist, args...)
end

function simulate(mvn::MvNormal, n_exp::Integer, horizon::Integer)
    dist = Distributions.MvNormal(mvn.mu, mvn.Sigma)
    d = length(mvn.mu)
    out = zeros(d, n_exp, horizon)
    for t=1:horizon
        out[:,:,t] = rand(dist,n_exp)
    end
    return out
end


discretize(dmp::DiscreteMarkovProcess) = dmp

# function discretize(dmp::DiscreteMarkovProcess)
#     nodes = dmp.values
#     n_nodes = size(nodes, 1)
#     integration_nodes = [nodes for i=1:n_nodes]
#     integration_weights = [dmp.transitions[i, :] for i=1:n_nodes]
#     return DiscretizedProcess(nodes, integration_nodes, integration_weights)
# end

#type DiscretizedContinousProcess <: AbstractDiscretizedProcess
#    nodes::Array{Float64,2}
#    integration_nodes::Array{Float64,3}
#    integration_weights::Array{Float64,1}
#end

type VAR1 <: ContinuousProcess
    M::Array{Float64,1}
    R::Array{Float64,2}
    Sigma::Array{Float64,2}
end

# Questions:
# Check :  Only for myVAR.M = 0 and myVAR.R and myVAR.Sigma is diagonal...?
# If that works add for M= 0

"""
    discretize(VAR_info::myVAR, n_states::Float64, n_integration::Float64, n_std::Int64=3, μ=::Int64=0)

Discretize a VAR(1)
This function takes:
  - VAR_info -  an object of the type "VAR1.myvar". It contains 3 arrays: VAR_info.M - mean; VAR_info.R - autocorrelation matrices; VAR_info.Sigma - covariance matrix of the innovations;
  - S - number of states for endogenous processes;
  - S_exo - number of states for exogenous processes;
  - n_std : int, optional(default=3). The number of standard deviations to each side the processes
   should s
The function returns:
  nodes::Matrix{Float64}
  integration_nodes::Array{Matrix{Float64}}
  integration_weights::Array{Vector{Float64}}
"""

function VAR1(R::Array{Float64,2}, Sigma::Array{Float64,2})
    M = zeros(size(R, 1))
    return VAR1(M, R, Sigma)
end

function discretize(var::VAR1)
    d = size(var.R, 1)
    n_states = ones(Int64, d)*5
    n_integration = ones(Int64, d)*5
    dis = discretize(var, n_states, n_integration)
    return dis
end


function discretize(var::VAR1, n_states::Array{Int64,1}, n_integration::Array{Int64,1}; n_std::Int64=2,)

    R = var.R
    M = var.M
    Sigma = var.Sigma
    S = QE.solve_discrete_lyapunov(R,Sigma)  # asymptotic variance
    sig = diag(S)
    min = var.M - n_std*(sig)
    max = var.M + n_std*(sig)
    grid = Dolo.CartesianGrid(min,max,n_states)
    # discretize innovations
    x,w = QE.qnwnorm(n_integration, zeros(size(var.Sigma,1)), var.Sigma)
    integration_nodes = [ cat(1,[(M + R*(Dolo.node(grid, i)-M) + x[j,:])' for j=1:size(x,1)]...) for i in 1:Dolo.n_nodes(grid)]
    integration_weights = [w for i in 1:Dolo.n_nodes(grid)]
    return DiscretizedProcess(grid, integration_nodes, integration_weights)
end
#
# function discretize(VAR_info::VAR1, n_states::Array{Int64,1}, n_integration::Array{Int64,1}; n_std::Int64=3, μ::Int64=0)
#
#       ###############################################################################
#
#       n  = size(VAR_info.M, 1);
#       # State space
#       # Collect all the state values for each variable in VAR
#       states = collect(zeros(n_states[ii]) for ii in 1:n);
#       a_bar=zeros(n)
#       # That is true iff R and Sigma is diagonal... better to use simulate?
#       a_bar = [n_std * sqrt(VAR_info.Sigma[j, j]^2 / (1 - VAR_info.R[j, j]^2)) for j in 1:n]
#       y = [linspace(-a_bar[j], a_bar[j], n_states[j]) for j in 1:n]
#       states = [collect(y[j]) for j in 1:n]
#       # Take cartesian product to have state space for endogenuous part
#       yⁱ= collect(Base.product([states[i] for i in 1:size(states, 1)]...))
#       nodes = cat( 1, [[e...]' for e in yⁱ]... )
#
#       x = collect(zeros(n_integration[ii]) for ii in 1:n);
#       w = collect(zeros(n_integration[ii]) for ii in 1:n);
#       for j in 1:n
#           x[j], w[j] = VAR_info.Sigma[j, j]*gauss(Float64, n_integration[j])
#       end
#       # Computing the nodes of the exo process
#       Eⁱʲ=collect(zeros(n_integration[ii]) for ii in 1:n);
#       #Compute states of exogenous process
#       # j stands for a variables in the VAR
#
#       [Eⁱʲ[j] =  1/(2*VAR_info.Sigma[j, j]*sqrt(pi)).*exp(-(x[j]- μ).^2./(2*VAR_info.Sigma[j, j]^2))  for j = 1:n]
#       #Compute Cartesian product
#
#       Eⁱʲ_c = collect(Base.product([Eⁱʲ[i] for i in 1:size(Eⁱʲ, 1)]...))
#       #Compute future states for endogenous variables
#       # kron([n_integration[i] for i in 1:n]...) computes the possible combination of states of innovations
#       if n>2
#         yⁱʲ  = zeros(kron([n_integration[i] for i in 1:n]...), n, kron([n_states[i] for i in 1:n]...));
#
#         for ii in 1:kron([n_states[i] for i in 1:n]...)
#
#             [ yⁱʲ[ss, :, ii] = VAR_info.M + VAR_info.R*collect(yⁱ[ii])+collect(Eⁱʲ_c[ss]) for ss = 1:kron([n_integration[i] for i in 1:n]...) ]
#
#         end
#         integration_nodes = collect(yⁱʲ[:, :, ii] for ii in 1:kron([n_states[i] for i in 1:n]...))
#
#         # Get transition probabilities
#
#         integration_weights = repmat(kron([w[i] for i in 1:n]...), 1, kron([n_integration[i] for i in 1:n]...))
#         integration_weights =collect(integration_weights[:, ii] for ii in 1:kron([n_integration[i] for i in 1:n]...))
#       else
#         yⁱʲ  = zeros(n_integration[n], n, n_states[n]);
#
#         for ii in 1:n_states[n]
#             [ yⁱʲ[ss, :, ii] = VAR_info.M + VAR_info.R*collect(yⁱ[ii])+collect(Eⁱʲ_c[ss]) for ss = 1:n_integration[n]]
#
#         end
#         integration_nodes = collect(yⁱʲ[:, :, ii] for ii in 1:n_states[n])
#
#         # Get transition probabilities
#
#         integration_weights = w[n]
#         integration_weights =[integration_weights]
#       end
# #  return nodes, integration_nodes, integration_weights
#     grid = CartesianGrid(-a_bar, a_bar, n_states)
#     return DiscretizedProcess(grid, integration_nodes[1], integration_weights[1])
# end



function simulate(var::VAR1, N::Int64, T::Int64)

  """
  This function takes:
    - VAR_info -  an object of the type "VAR1.myvar"
    which contains 3 arrays: mean and autocorrelation matrices of the VAR(1) process
      and covariance matrix of the innovations;
    - N - number of simulations.
    - T - number of periods to simulate;

  The function returns:
    - XN - simulated data, ordered as [n, N, T], where n is the number of variables;
    #- E - simulated innovations, ordered as [N, T, n], where n is the number of variables;
  """
  n = size(var.M, 1);
  # srand(123) # Setting the seed
  d = MvNormal(var.Sigma);
  XN=zeros(N, T, n);
  E=zeros(N, T, n);
  #Inicial conditions
  XN[:, 1, :]=rand(N, n);
  for jj = 1:1:N
        E[jj, :, :] = rand(d, T)';
        for ii =1:1:T-1
            XN[jj, ii+1, :]=var.M+var.R*(XN[jj, ii, :]-var.M)+E[jj, ii, :];
        end
  end
  return permutedims(XN,[3,1,2])  #, E
end





function ErgodDist(VAR_info::VAR1, T::Int64, N::Int64)

       # Simulate for T period, N times;
       #  #    - simulate a VAR1 process for T periods N times#
       n = size(VAR_info.M, 1);
       srand(123) # Setting the seed
       d = MvNormal(VAR_info.Sigma);
       VAR_process=sim_VAR1.simulate_var(VAR_info::VAR1.myvar, T::Int64, N::Int64)
       E = VAR_process[2];
       XN = VAR_process[1];
       # Computing moments
       Mean_sim = mean(squeeze(mean(VAR_process[1], 2), 2), 1);
       #Computing the std of simulted Processes (Covariance matrix)
       E1_d = (E[:, :, 1] - repeat(mean(E[:, :, 1], 2), inner=[1, T]));
       cov(E1_d[:, 1])    # which is X1_d[:, 1]'*X1_d[:, 1]/T
       diag(cov(E1_d[:, :]))   # covariances across simulation
       Sigma_sim = mean(diag(cov(E1_d[:, :])) )  # mean of covariances across simulations
       # Autocorrelation matrix
       X_d=zeros(N, T, 2)
       X_d0=zeros(N, T-1, 2);
       X_d1=zeros(N, T-1, 2);
       R_sim = zeros(2, 2);
       for ii in [1, 2]
       X_d[:, :, ii] = (XN[:, :, 2] - repeat(mean(XN[:, :, ii], 2), inner=[1, T]));
       X_d0[:, :, ii]  = X_d[:, 1:end-1, ii];
       X_d1[:, :, ii] = (XN[:, 2:end, ii] - repeat(mean(XN[:, 2:end, ii], 2), inner=[1, T-1]));
       R_sim[ii, ii] =  mean(diag( cov(X_d0[:, :, ii], X_d1[:, :, ii] )/sqrt(var(X_d0[:, :, ii]))/sqrt(var(X_d1[:, :, ii]))  ))
         if ii == 2
             R_sim[ii-1, ii]  =   mean(diag( cov(X_d0[:, :, ii-1], X_d1[:, :, ii] )/sqrt(var(X_d0[:, :, ii-1]))/sqrt(var(X_d1[:, :, ii])) ))
             R_sim[ii, ii-1]  =   mean(diag( cov(X_d0[:, :, ii], X_d1[:, :, ii-1] )/sqrt(var(X_d0[:, :, ii]))/sqrt(var(X_d1[:, :, ii-1])) ))
         end
       end

       return Mean_sim, Sigma_sim, R_sim
end

# compatibility names
typealias AR1 VAR1
typealias MarkovChain DiscreteMarkovProcess
typealias Normal MvNormal
