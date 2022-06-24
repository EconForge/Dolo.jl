using Dolo


model = Dolo.yaml_import("examples/models/consumption_savings_iid.yaml")


res =Dolo.get_factory(model,"expectation")

sol_direct = Dolo.time_iteration_direct(model, maxit=500;tol_η=1e-12)
sol_iti = Dolo.improved_time_iteration(model, maxit=500; dr0=sol_direct.dr)
sol = Dolo.time_iteration(model; dr0=sol_direct.dr)

tab_d = Dolo.tabulate(model, sol_direct.dr, :w)
tab_i = Dolo.tabulate(model, sol_iti.dr, :w)
tab = Dolo.tabulate(model, sol.dr, :w)

pl = plot(tab_d[:w], tab_d[:c])
plot!(tab[:w], tab_i[:c])
plot!(tab[:w], tab[:c])
# plot!(tab[:w],tab[:w])
pl
