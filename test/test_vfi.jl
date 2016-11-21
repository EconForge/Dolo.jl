# d = Pkg.dir("Dolo")
# cd(d)

    using Dolo

    model = yaml_import("examples/models/rbc_dtcc_mc.yaml")
    @time dr = Dolo.time_iteration(model, verbose=true, maxit=2000)


model

val0 = Dolo.evaluate_policy(model, dr, maxit=5000)

kvec = collect(linspace(5,30,10))
vec0 = [Dolo.evaluate(val0, 1, [k])[1] for k in kvec]

n_vec = [Dolo.evaluate(dr, 1, [k])[1] for k in kvec]
i_vec = [Dolo.evaluate(dr, 1, [k])[2] for k in kvec]
v_vec = [Dolo.evaluate(dr, 1, [k])[3] for k in kvec]

v_vec


n_vec, i_vec

include(joinpath("..","src","vfi.jl"))
import temp

drv, controls = temp.solve_policy(model, dr)


initial_x, lower, upper, drv = temp.solve_policy(model, dr)
vec = [Dolo.evaluate(drv, 1, [k])[1] for k in kvec]


using DataFrames
methods(DataFrame)
# df = DataFrame(Dict(:k=>kvec,:v1=>vec0,:v2=>vec))
plot(df,x=:k,y=:v1,y2=:v2)
using Gadfly
l1 = layer(x=kvec,y=vec0,Geom.line)
l2 = layer(x=kvec,y=vec,Geom.line)
plot(l1,l2)




kvec
vec
vec0



upper

s = Dolo.nodes(drv.grid)[25,:]
x0 = controls[1][25,:]

beta = 0.96

dprocess = Dolo.discretize(model.exogenous)
i = 1
p = model.calibration[:parameters]
x0

evaluate(drv,i,x0)

ivec = linspace(x0[2]*0.5,x0[2]*3)
xvec = [[x0[1],xx] for xx in ivec]
resp_vec = [temp.update_value(model, beta, dprocess, drv, i, s, xx, p) for xx in xvec]

evaluate(drv,i,x0)

l1 = layer(x=ivec,y=resp_vec,Geom.line)
l2 = layer(xintercept=[x0[2]], Geom.vline)
l3 = layer(yintercept=[evaluate(drv,i,x0)[1]], Geom.hline)

plot(l1,l2,l3)

x0

vec0 - vec0


using Gadfly
kvec
l1 = layer(x=kvec, y=vec, Geom.line)
l0 = layer(x=kvec, y=vec0, Geom.line)
plot(l0,l1)

plot(x=kvec, y=vec)

vec
vec0
=|?"""""""""""""""""""""""""""""""+_"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""l1 = layer(x=kvec, y=vec)[-+03
.]
