
using Dolo

fname = "examples/rbc_dtmscc.yaml"

model = yaml_import(fname)

residuals(model) # check that residuals are indeed correct

# check some calibrated values
qeval(s) = eval_with(model.calibration, s)
qeval("k")
qeval("k*0.5")
qeval("k*1.5")


# solution

# initial guess is equal to steady-state investment/labour
phi_0 = constant_guess(model)

s = model.calibration[:states]
evaluate(phi_0, 1, repmat(s', 5,1))


# value function iteration will need an initial guess for the value:
@time initval = evaluate_policy(model, phi_0)

# note that fmincon seems to fail when started very close to optimum
# for this reason we limit the number of iterations
@time controls, value = solve_policy(model, phi_0, initval, Dict(:maxit=>200))


# one can also iterate on euler equations
@time time_iteration(model, phi_0)



using Gadfly
x = collect(linspace(qeval("k*0.5"), qeval("k*1.5"),size(controls,2)))
Gadfly.plot(x=x, y=slice(value,1,:))


l1a = layer(x=x,y=controls[1,:,1], Geom.line)
l1b = layer(x=x,y=controls[2,:,1], Geom.line)
l2a = layer(x=x,y=controls[1,:,2], Geom.line)
l2b = layer(x=x,y=controls[2,:,2], Geom.line)
lx = layer(x=x,y=x*qeval("delta"), Geom.line)
p = layer(x=[qeval("k")],y=[qeval("delta*k")],Geom.point)
plot(l1a,l1b,l2a,l2b,lx,p)
