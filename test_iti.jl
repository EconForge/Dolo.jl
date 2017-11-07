import Dolo

import Dolo: n_nodes, n_inodes, nodes, CachedDecisionRule
import Dolo: invert_jac

model = Dolo.yaml_import("examples/models/rbc_dtcc_mc.yaml")
# model = Dolo.yaml_import("examples/models/sudden_stop.yaml")
# dri = Dolo.time_iteration(model)
# Dolo.nodes( model.grid )

@time dr = Dolo.improved_time_iteration(model, verbose=false)

@time dr = Dolo.improved_time_iteration(model, verbose=false, method=:gmres)

@time dr = Dolo.improved_time_iteration(model, verbose=false, method=:iti)

R_i, D_i, J_ij_0, S_ij, dumdr = Dolo.improved_time_iteration(model)
J_ij = J_ij_0*(-1)

maxabs(R_i)
rsol, it, lam, errors = invert_jac(R_i, D_i, J_ij, S_ij, dumdr)

π_i, M_ij, S_ij = Dolo.preinvert(R_i, D_i, J_ij, S_ij)

x0 = deepcopy(π_i)



Π_i = deepcopy(π_i)
Dolo.d_filt_dx!(Π_i, π_i, M_ij, S_ij, dumdr)

L = LinearThing(M_ij, S_ij, dumdr)

import Dolo: ListOfPoints



function app(L::LinearThing,x::Vector{ListOfPoints{n_x}}) where n_x
   xx = deepcopy(x)
   Dolo.d_filt_dx!(xx, x, L.M_ij, L.S_ij, L.I)
   return xx
end

@time xx = app(L,x0)


maxabs(xx)
@time xx = app(L,xx)





function solveit(L::LinearThing,v::AbstractVector{Float64},tol=1e-8)
   x = copy(v)
   tot = copy(x)
   maxit= 1000
   it=0
   while it<maxit
      it+=1
      x = mul(L,x)
      tot += x
   end
   return tot
end

size(L,1)
v = rand(size(L,1))

w = solveit(L,v)

maximum(abs,L*w-v)


@time rsol, it, lam, errors = invert_jac(R_i, D_i, J_ij, S_ij, dumdr)




n_m = length(x0)
N = length(x0[1])
n_x = length(x0[1][1])
dd = rand(n_x,N,n_m)
v = rand(n_x,N,n_m)[:]
w =rand(n_x,N,n_m)[:]

x0

r0 = L*x0
rr = L*dd

@time res = L*v

@time sol = gmres(L, v, verbose=true, tol=1e-10, maxiter=1000, restart=100)

sol0 = solveit(L, v)
@time rsol, it, lam, errors = invert_jac(R_i, D_i, J_ij, S_ij, dumdr)



@time sol = gmres(L, v)

sol - sol0





@time sol2 = solveit(L, v)

sol - sol2



L*sol - v





L*v


res = L*dd

@time size(res)

dprocess = Dolo.discretize(model.exogenous)
grid = model.grid
init_dr = Dolo.ConstantDecisionRule(model.calibration[:controls])

parms = model.calibration[:parameters]

#  n_m = Dolo.n_nodes(dprocess) # number of exo states today
n_m = max(n_nodes(dprocess), 1) # number of exo states today
n_mt = n_inodes(dprocess,1)  # number of exo states tomorrow
n_s = length(model.symbols[:states]) # number of endo states

s = nodes(grid)
N_s = size(s,1)
n_x = size(model.calibration[:controls],1)
#  N_m = Dolo.n_nodes(dprocess) # number of grid points for exo_vars

#  x0 = [repmat(model.calibration[:controls]',N_s) for i in 1:N_m] #n_x N_s n_m
x0 = [init_dr(i, s) for i=1:n_m]
ddr=CachedDecisionRule(dprocess, grid, x0)
ddr_filt =
L*sol - v





L*v


res = L*dd

@time size(res)

dprocess = Dolo.discretize(model.exogenous)
grid = model.grid
init_dr = Dolo.ConstantDecisionRule(model.calibration[:controls])

parms = model.calibration[:parameters]

#  n_m = Dolo.n_nodes(dprocess) # number of exo states today
n_m = max(n_nodes(dprocess), 1) # number of exo states today
n_mt = n_inodes(dprocess,1)  # number of exo states tomorrow
n_s = length(model.symbols[:states]) # number of endo states

s = nodes(grid)
N_s = size(s,1)
n_x = size(model.calibration[:controls],1)
#  N_m
L*sol - v





L*v


res = L*dd

@time size(res)

L*sol - v





L*v


res = L*dd

@time size(res)

dprocess = Dolo.discretize(model.exogenous)
grid = model.grid
init_dr = Dolo.ConstantDecisionRule(model.calibration[:controls])

parms = model.calibration[:parameters]

#  n_m = Dolo.n_nodes(dprocess) # number of exo states today
n_m = max(n_nodes(dprocess), 1) # number of exo states today
n_mt = n_inodes(dprocess,1)  # number of exo states tomorrow
n_s = length(model.symbols[:states]) # number of endo states

s = nodes(grid)
N_s = size(s,1)
n_x = size(model.calibration[:controls],1)
#  N_m = Dolo.n_nodes(dprocess) # number of grid points for exo_vars

#  x0 = [repmat(model.calibration[:controls]',N_s) for i in 1:N_m] #n_x N_s n_m
x0 = [init_dr(i, s) for i=1:n_m]
ddr=CachedDecisionRule(dprocess, grid, x0)
ddr_filt = CachedDecisionRule(dprocess, grid, x0)
Dolo.set_va
dprocess = Dolo.discretize(model.exogenous)
grid = model.grid
init_dr = Dolo.ConstantDecisionRule(model.calibration[:controls])

parms = model.calibration[:parameters]

#  n_m = Dolo.n_nodes(dprocess) # number of exo states today
n_m = max(n_nodes(dprocess), 1) # number of exo states today
n_mt = n_inodes(dprocess,1)  # number of exo states tomorrow
n_s = length(model.symbols[:states]) # number of endo states

s = nodes(grid)
N_s = size(s,1)
n_x = size(model.calibration[:controls],1)
#  N_m = Dolo.n_nodes(dprocess) # number of grid points for exo_vars

#  x0 = [repmat(model.calibration[:controls]',N_s) for i in 1:N_m] #n_x N_s n_m
x0 = [init_dr(i, s) for i=1:n_m]
ddr=CachedDecisionRule(dprocess, grid, x0)
ddr_filt = CachedDecisionRule(dprocess, grid, x0)
Dolo.set_va= Dolo.n_nodes(dprocess) # number of grid points for exo_vars

#  x0 = [repmat(model.calibration[:controls]',N_s) for i in 1:N_m] #n_x N_s n_m
x0 = [init_dr(i, s) for i=1:n_m]
ddr=CachedDecisionRule(dprocess, grid, x0)
ddr_filt = CachedDecisionRule(dprocess, grid, x0)
Dolo.set_vaCachedDecisionRule(dprocess, grid, x0)
Dolo.set_values!(ddr,x0)

import Dolo: to_LOP
using StaticArrays

s_ = to_LOP(s)
x_ = [to_LOP(el) for el in x0]
p_ = SVector(parms...)

O(u) = Dolo.euler_residuals(model,s_,u,ddr,dprocess,p_,with_jres=false)

epsilon=1e-8
x0 = copy(x_)
R0 = O(x0)

import Dolo: ListOfPoints, to_LOJ


J[1][1]


# @time sol = Dolo.time_iteration(model, maxit=1000)


@time ret = Dolo.improved_time_iteration(model, complementarities=true)


# @time ddr = Dolo.improved_time_iteration(model, complementarities=true)

vvvv = Dolo.CachedDecisionRule{Dolo.CubicDR{Dolo.EmptyGrid,Dolo.CartesianGrid{2},2,2},Dolo.DiscretizedIIDProcess}




(::Int64, ::Int64, ::Array{SVector{2,Float64},1})
