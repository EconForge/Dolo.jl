import Dolo

model = Dolo.Model(joinpath(Dolo.pkg_path, "examples","models","rbc_mc.yaml"))

Dolo.time_iteration(model)
dr0 = Dolo.time_iteration(model, tol_η=1e-6; grid=Dict(:n=>[5])).dr
sol = Dolo.time_iteration(model, maxit=1000; grid=Dict(:n=>[100]))

# faster version: limit the number of steps in the newton solver
@time sol = Dolo.time_iteration(model, maxit=1000; solver=Dict(:maxit=>2), grid=Dict(:n=>[100]))

# or use explicit formulas given in the model
@time sol = Dolo.time_iteration(model, maxit=1000; solver=Dict(:type=>:direct), grid=Dict(:n=>[100]))
