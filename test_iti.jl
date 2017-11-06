import Dolo

import Dolo: n_nodes, n_inodes, nodes, CachedDecisionRule
import Dolo: invert_jac

model = Dolo.yaml_import("examples/models/rbc_dtcc_iid.yaml")
dri = Dolo.time_iteration(model)
Dolo.nodes( model.grid )
R_i, D_i, J_ij, S_ij, dumdr = Dolo.improved_time_iteration(model)

R_i


inv( D_i[1][90] )*D_i[1][90]

rsol, it, lam, errors = invert_jac(R_i, D_i, J_ij, S_ij, dumdr)

errors

D_i[1][50]\R_i[1][50]
D_i



π_i, M_ij, S_ij = Dolo.preinvert(R_i, D_i, J_ij, S_ij)

x0 = deepcopy(π_i)


struct LinearThing
   M_ij
   S_ij
   I
end

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


@time xx = app(L,xx)
maxabs(xx)

L*x0

@time

import Base.size
import Base.eltype
import Base.*
using StaticArrays

eltype(L::LinearThing) = Float64
function shape(L::LinearThing)
   n_m = size(L.M_ij,1)
   N = size(L.M_ij[1,1],1)
   n_x = size(L.M_ij[1,1][1],1)
   return (n_x, N, n_m)
end
size(L::LinearThing,d) = prod(shape(L))



function *(L::LinearThing,x::Vector{ListOfPoints{n_x}}) where n_x
   xx = deepcopy(x)
   Dolo.d_filt_dx!(xx, x, L.M_ij, L.S_ij, L.I)
   return xx
end

function *(L::LinearThing,m::Array{Float64, 3})
   n_x,N,n_m = size(m)
   x = [reinterpret(SVector{n_x,Float64},m[:,:,i],(N,)) for i=1:n_m]
   y = deepcopy(x)
   xx = L*y
   rr = [reinterpret(Float64, xx[i], (n_x,N)) for i=1:length(xx)]
   rrr = cat(3,rr...)
   return reshape(rrr, n_x,N,n_m)
end

function *(L::LinearThing,v::AbstractVector{Float64})
   m = copy(v)
   sh = shape(L)
   mm = reshape(m, sh...)
   mmm = L*mm + mm
   return mmm[:]
end

function mul(L::LinearThing,v::AbstractVector{Float64})
   m = copy(v)
   sh = shape(L)
   mm = reshape(m, sh...)
   mmm = L*mm
   return mmm[:]
end


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

mul(L,v)

function powww(L::LinearThing, v)
   x = copy(v)
   lams = []
   err0 = 1.0
   for i=1:100
      x = mul(L,x)
      err = maximum(abs(x))
      lam = err/err0
      push!(lams, lam)
   end
   return lams
end

ll = powww(L,v)

ll


powm(L)

@time tot = solveit(L,v)

mul(L,tot)-v





import Base.A_mul_B!

function A_mul_B!(w::AbstractVector{Float64},L::LinearThing,v::AbstractVector{Float64})
   w[:] = L*v
   return w
end

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

import IterativeSolvers

@time sol = gmres(L, v, verbose=true, tol=1e-10, maxiter=1000, restart=100)







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
