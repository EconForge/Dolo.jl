using PyCall

@pyimport dolo

modelpy = dolo.yaml_import("examples/models/consumption_savings_iid.yaml")
modelpy.__exogenous__ = modelpy.exogenous.processes[1]

@pyimport dolo.algos.time_iteration_helpers as tih
Fpy = tih.Euler(modelpy)


using Dolo

model = Dolo.yaml_import("examples/models/consumption_savings_iid.yaml")


Dolo.time_iteration(model; maxit=500)



F = Dolo.Euler(model)

# set equal dprocess

Dolo.n_inodes(F.dprocess,1)
Dolo.inode(F.dprocess,1,1)
Dolo.iweight(F.dprocess,1,1)

Fpy.dprocess.n_inodes(1)
Fpy.dprocess.iweight(1,0)

jpoints = [Dolo.inode(F.dprocess, 1, j) for j=1:Dolo.n_inodes(F.dprocess, 1)]
jweights = [Dolo.iweight(F.dprocess, 1, j) for j=1:Dolo.n_inodes(F.dprocess, 1)]

ppoints  = [Fpy.dprocess.inode(0,j) for j=0:6]
pweights  = [Fpy.dprocess.iweight(0,j) for j=0:6]

F.dprocess.integration_weights = pweights
F.dprocess.integration_nodes = copy(hcat(ppoints...)')



# check derivatives without complementarities

flatten(x0::Dolo.MSM) = cat([e for e in x0.data]...;dims=1)


# same residual
x0_py = Fpy.x0
respy = Fpy(x0_py, x0_py)
resj = F(F.x0, F.x0; ignore_constraints=false)
maximum(abs, (respy.data)[:] - flatten(resj))


# same jacobian
rpy, Jpy = Fpy.d_A(x0_py, x0_py)
Jj = Dolo.df_A(F, F.x0, F.x0)
diff = [Jj.data[n] - (Jpy.data)[n,:,:] for n=1:length(Jj.data)]
ddiff = cat(diff...; dims=1)
maximum(abs, ddiff)


# same step
delta_py = Jpy.solve(rpy)
delta_jl = Jj\resj
maximum(abs, (delta_py.data)[:] - flatten(delta_jl))

x0j = F.x0
x0p = Fpy.x0


x0j = x1j
x0p = x1p

x1j = x0j - delta_jl
x1p = x0p - delta_py
maximum(abs, (x1p.data)[:] - flatten(x1j))


res1j = F(x1j, x0j; ignore_constraints=false)
res1p = Fpy(x1p, x0p; ignore_constraints=false)

J1j = Dolo.df_A(F, x1j, x0j)
r1p, J1p = Fpy.d_A(x1p, x0p)

Dolo.norm(res1j)
res1p.norm()

maximum(abs, (res1p.data)[:] - flatten(res1j))


diff1 = [J1j.data[n] - (J1p.data)[n,:,:] for n=1:length(J1j.data)]
ddiff1 = cat(diff...; dims=1)
maximum(abs, ddiff1)

#####



time_iteration(model, maxit=100)

dolo.time_iteration(modelpy, maxit=100)

sol = Dolo.improved_time_iteration(model)


tab = Dolo.tabulate(model, sol.dr, :w)

using SimplePlots

plot(tab[:w], tab[:w], ylims=(0,2))
plot!(tab[:w], tab[:c])

