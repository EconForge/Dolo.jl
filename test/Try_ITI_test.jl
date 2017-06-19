
path = Pkg.dir("Dolo")

# Pkg.build("QuantEcon")
import Dolo
import Bruteforce_module

###############################################################################
## Markov Chain
###############################################################################
filename = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
# model = Dolo.Model(Pkg.dir("Dolo", "examples", "models", "rbc_dtcc_mc.yaml"), print_code=true)
model = Dolo.yaml_import(filename)

@time dr_ITI  = Bruteforce_module.improved_time_iteration(model; verbose=true, tol = 1e-06, smaxit=50)
# @time dr_ITI_2  = Bruteforce_module.improved_time_iteration(model, dprocess,init_dr)

dprocess = Dolo.discretize( model.exogenous )
init_dr = Dolo.ConstantDecisionRule(model.calibration[:controls])
maxbsteps= 50
complementarities = false


parms = model.calibration[:parameters]

#  n_m = Dolo.n_nodes(dprocess) # number of exo states today
n_m = max(Dolo.n_nodes(dprocess), 1) # number of exo states today
n_mt = Dolo.n_inodes(dprocess,1)  # number of exo states tomorrow
n_s = length(model.symbols[:states]) # number of endo states

grid = Dolo.Dict()
grid = Dolo.get_grid(model, options=grid)
s = Dolo.nodes(grid)
N_s = size(s,1)
n_x = size(model.calibration[:controls],1)
#  N_m = Dolo.n_nodes(dprocess) # number of grid points for exo_vars

#  x0 = [repmat(model.calibration[:controls]',N_s) for i in 1:N_m] #n_x N_s n_m
x0 = [init_dr(i, s) for i=1:n_m]
ddr=Dolo.CachedDecisionRule(dprocess, grid, x0)
ddr_filt = Dolo.CachedDecisionRule(dprocess, grid, x0)
Dolo.set_values!(ddr,x0)

steps = 0.5.^collect(0:maxbsteps)

if complementarities == true
 x_lb = Array{Float64,2}[cat(1, [Dolo.controls_lb(model, Dolo.node(dprocess, i), s[n, :], parms)' for n=1:N_s]...) for i=1:n_m]
 x_ub = Array{Float64,2}[cat(1, [Dolo.controls_ub(model, Dolo.node(dprocess, i), s[n, :], parms)' for n=1:N_s]...) for i=1:n_m]
end


x=x0
## memory allocation
jres = zeros(n_m,n_mt,N_s,n_x,n_x)
S_ij = zeros(n_m,n_mt,N_s,n_s)

######### Loop     for it in range(maxit):
it=0
it_invert=0

# res_init = euler_residuals(model,s,x,ddr,dprocess,parms ,set_dr=false, jres=jres, S_ij=S_ij)

x

N_s = size(s,1) # Number of gris points for endo_var
n_s = size(s,2) # Number of states
n_x = size(x[1],2) # Number of controls

# P = dprocess.values
# Q = dprocess.transitions

# n_ms = Dolo.n_nodes(dprocess)  # number of exo states today
n_ms = size(x,1)  # number of exo states today
n_mst = Dolo.n_inodes(dprocess,1)  # number of exo states tomorrow
# n_mv = size(Dolo.node(dprocess, 1),1)  # number of exo variable

res = zeros(n_ms, N_s, n_x)
with_jres = false
if with_jres == true
  if jres== nothing
    jres = zeros((n_ms,n_mst,N_s,n_x,n_x))
  end
  if S_ij== nothing
    S_ij = zeros((n_ms,n_mst,N_s,n_s))
  end
end



Dolo.node(dprocess,i_ms)
[Dolo.node(dprocess,i_ms) for i in 1:N_s]
repmat(Dolo.node(dprocess,i_ms),1)
i_ms=1
m_prep = [repmat(Dolo.node(dprocess,i_ms),1) for i in 1:N_s]
m=(hcat([e' for e in m_prep]...))'
I_ms =1
M_prep = [repmat(Dolo.inode(dprocess, i_ms, I_ms), 1) for i in 1:N_s]
M=(hcat([e' for e in M_prep]...))'
w = Dolo.iweight(dprocess, i_ms, I_ms)

[s for i in 1:2]
repmat(s,2)
s = Dolo.nodes(grid)
S = Dolo.transition(model, m, repmat(s,2), x[i_ms], M, parms)




S = Dolo.transition(model, m, s, vec(x[i_ms]), M, parms)

vec(x[i_ms])






t
