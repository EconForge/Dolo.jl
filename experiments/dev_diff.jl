using Dolo

model = Model("examples/models/neoclassical.yaml")

sol = improved_time_iteration(model)


using StaticArrays

m, s, x , p = model.calibration[:exogenous, :states, :controls, :parameters]
m = SVector(m...)
s = SVector(s...)
x = SVector(x...)
p = SVector(p...)

using ForwardDiff

sol.dr(1,s)

a = sol.dr.grid_endo.min
b = sol.dr.grid_endo.max
n = sol.dr.grid_endo.n

f = u->Dolo.splines.eval_UC_spline(a,b,n,sol.dr.itp[1],u)


Dolo.splines.eval_UC_spline_(a,b,n,sol.dr.itp[1],s)

f(s)

ForwardDiff.jacobian(f, s)

sol.dr(1,s)


ForwardDiff.jacobian(u->sol.dr(1,u), s)

V,dV = sol.dr(Val((0,2)), 1, [s, s, s, s])

f = u->arbitrage(model, m, s, x, m, s, u, p)

f(x)

autodiff = ForwardDiff.jacobian(f,x)

symbolic = arbitrage(model, Val(6), m, s, x, m, s, x, p)

autodiff - symbolic


F = Dolo.Euler(model)

J = Dolo.df_A(F, F.x0, F.x0)




@time sol = improved_time_iteration(model)
