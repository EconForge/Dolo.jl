import QuantEcon
const QE = QuantEcon
import Base:rand
import Distributions

abstract AbstractExogenous

abstract AbstractProcess <: AbstractExogenous
abstract DiscreteProcess <: AbstractProcess
abstract ContinuousProcess <: AbstractProcess
abstract IIDExogenous <: AbstractProcess

abstract AbstractDiscretizedProcess

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
node(dp::DiscretizedProcess, i) = node(dp.grid,i)
n_inodes(dp::DiscretizedProcess, i::Int) = size(dp.integration_nodes[i], 1)
inode(dp::DiscretizedProcess, i::Int, j::Int) = dp.integration_nodes[i][j, :]
iweight(dp::DiscretizedProcess, i::Int, j::Int) = dp.integration_weights[i][j]

# date-t grid is unstructured
type DiscreteMarkovProcess <: AbstractDiscretizedProcess
    grid::UnstructuredGrid
    transitions::Matrix{Float64}
    values::Matrix{Float64}
end

DiscreteMarkovProcess(transitions::Matrix{Float64}, values::Matrix{Float64}) =
    DiscreteMarkovProcess(UnstructuredGrid(values), transitions, values)

n_nodes(dp::DiscreteMarkovProcess) = size(dp.values, 1)
n_inodes(dp::DiscreteMarkovProcess, i::Int) = size(dp.values, 1)
inode(dp::DiscreteMarkovProcess, i::Int, j::Int) = dp.values[j, :]
iweight(dp::DiscreteMarkovProcess, i::Int, j::Int) = dp.transitions[i,j]
node(dp::DiscreteMarkovProcess, i) = dp.values[i, :]


# date-t grid is empty

type DiscretizedIIDProcess <: AbstractDiscretizedProcess
    # nodes::Matrix{Float64}
    grid::EmptyGrid
    integration_nodes::Matrix{Float64}
    integration_weights::Vector{Float64}
end

DiscretizedIIDProcess(x,w) = DiscretizedIIDProcess(EmptyGrid(), x, w)

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

function simulate(mvn::MvNormal, N::Integer, T::Integer, e0::Vector{Float64}; stochastic=true)
    dist = Distributions.MvNormal(mvn.mu, mvn.Sigma)
    d = length(mvn.mu)
    out = zeros(d, N, T)
    for i=1:N
        out[:,i,1] = e0
    end
    if stochastic
        for t =2:T
            out[:,:,t] = rand(dist, N)
        end
    end
    return out
end

function simulate(mvn::MvNormal, N::Integer, T::Integer; stochastic=true)
    e0 = zeros(size(mvn.mu,1))
    return simulate(mvn, N, T, e0; stochastic=stochastic )
end

function response(mvn::MvNormal, T::Integer, e1::AbstractVector)
    d = length(mvn.mu)
    out = zeros(d,T)
    out[:,2] = e1
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



function VAR1(R::Array{Float64,2}, Sigma::Array{Float64,2})
    M = zeros(size(R, 1))
    return VAR1(M, R, Sigma)
end

function discretize(var::VAR1)
    d = size(var.R, 1)
    n_states = ones(Int, d)*5
    n_integration = ones(Int, d)*5
    dis = discretize(var, n_states, n_integration)
    return dis
end

function discretize(var::VAR1, n_states::Array{Int,1}, n_integration::Array{Int,1}; n_std::Int=2,)
    R = var.R
    M = var.M
    Sigma = var.Sigma
    S = QE.solve_discrete_lyapunov(R,Sigma)  # asymptotic variance
    sig = diag(S)
    min = var.M - n_std*sqrt(sig)
    max = var.M + n_std*sqrt(sig)
    grid = Dolo.CartesianGrid(min,max,n_states)
    # discretize innovations
    x,w = QE.qnwnorm(n_integration, zeros(size(var.Sigma,1)), var.Sigma)
    integration_nodes = [ cat(1,[(M + R*(Dolo.node(grid, i)-M) + x[j,:])' for j=1:size(x,1)]...) for i in 1:Dolo.n_nodes(grid)]
    integration_weights = [w for i in 1:Dolo.n_nodes(grid)]
    return DiscretizedProcess(grid, integration_nodes, integration_weights)
end



function simulate(var::VAR1, N::Int, T::Int, e0::Vector{Float64}; stochastic=true)

  """
  This function takes:
    - var -  a VAR1 process
    - N - number of simulations.
    - T - number of periods to simulate;

  The function returns:
    - XN - simulated data, ordered as (n, N, T), where n is the number of variables;
  """
  n = size(var.M, 1);
  # srand(123) # Setting the seed
  d = MvNormal(var.Sigma);
  XN=zeros(N, T, n);
  E=zeros(N, T, n);
  if stochastic
      for jj = 1:1:N
            E[jj, :, :] = rand(d, T)';
      end
  end
  # Initial conditions
  for i=1:N
      XN[i, 1, :]=e0
      for jj = 1:1:N
        for ii =1:1:T-1
            XN[jj, ii+1, :]=var.M+var.R*(XN[jj, ii, :]-var.M)+E[jj, ii, :];
        end
      end
  end
  return permutedims(XN,[3,1,2])  #, E
end


function simulate(var::VAR1, N::Int, T::Int; stochastic=true)
    return simulate(var, N, T, var.M; stochastic=stochastic)
end



function ErgodDist(VAR_info::VAR1, T::Int, N::Int)

       # Simulate for T period, N times;
       #  #    - simulate a VAR1 process for T periods N times#
       n = size(VAR_info.M, 1);
       srand(123) # Setting the seed
       d = MvNormal(VAR_info.Sigma);
       VAR_process=sim_VAR1.simulate_var(VAR_info::VAR1.myvar, T::Int, N::Int)
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
