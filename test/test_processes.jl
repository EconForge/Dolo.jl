import Dolo
import Distributions


# test IID-Normal
mu = [0.0, 0.0]
Sigma = [1.0 0.0; 0.0 2.0]
d = Dolo.MvNormal(mu, Sigma)
sim = Dolo.simulate(d, 1, 200, [0.1, -0.01])
sim = Dolo.simulate(d, 1, 200, [0.1, -0.01]; stochastic=false)
dp = Dolo.discretize(d)


# test VAR
var = Dolo.VAR1([0.0,0.0],[0.99 0.0; 0.0 0.08],eye(2)*0.01)


sim = Dolo.simulate(var, 5, 200)
sim = Dolo.simulate(var, 5, 200, [0.8, 0.8])
sim = Dolo.simulate(var, 1, 200, [0.8, 0.8]; stochastic=false)

dp = Dolo.discretize(var)


import PyPlot







import temp


# markov chain
values = [ 0   0;
           0.1 0;
           0.1 0.1 ]

transitions = [ 0.1 0.4 0.5;
                0.1 0.8 0.1;
                0.2 0.0 0.8 ]

mc = DiscreteMarkovProcess(transitions, values)
dproc_mc = discretize(mc)
dproc_mc.nodes
dproc_mc.integration_weights
dproc_mc.integration_weights
@assert n_nodes(dproc_mc) == 3
dproc_mc

# normal law
Sigma = [0.05 0.01; 0.01 0.1]
mvn = MvNormal(Sigma)
dproc_mvn = discretize(mvn)
dproc_mvn.integration_weights
dproc_mvn.integration_nodes
dproc_mvn
@assert n_nodes(dproc_mvn) == 0


###############################################################################
M = rand(2)
Sigma = [0.5 0.5; 0.5 0.5]
R = [0.9 0.1; 0 .8];
n_states= [5,10]
n_states::Array{Int,1}
n = 2
n_integration=[7,4]
# 4/ Create a VAR1 type (in the module). Add new functions that operate on this type. Add them to the test file.

VAR_info = VAR1(M,R,Sigma)
using Base.product
D_Var = discretize(VAR_info, n_states, n_integration)
D_Var.nodes
D_Var.integration_nodes
D_Var.integration_weights[2]
[sum(D_Var.integration_weights[j]) for j in 1:kron([n_integration[j] for j in 1:n]...)]


#

M = rand(1)
Sigma = [0.5]'
R = [0.9]'
n_states= [5]
n_states::Array{Int,1}
n = 2
n_integration=[7]
# 4/ Create a VAR1 type (in the module). Add new functions that operate on this type. Add them to the test file.

VAR_info = VAR1(M,R,Sigma)
using Base.product
D_Var = discretize(VAR_info, n_states, n_integration)
D_Var.nodes
D_Var.integration_nodes
D_Var.integration_weights[2]
[sum(D_Var.integration_weights[j]) for j in 1:kron([n_integration[j] for j in 1:n]...)]
