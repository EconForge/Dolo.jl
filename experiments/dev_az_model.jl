using Dolo
using StaticArrays
using Plots

model = yaml_import("examples/models/az_model.yaml")

gr = Dolo.discretize(model)

# dr_pert = perturb(model)

dr_global = time_iteration(model,verbose=true)

sim = simulate(model, dr_global.dr, T=450, s0=[0.5])

p4 = plot(sim[1,:k,:],sim[1,:i,:], xlabel = "k", ylabel = "i")

p1 = plot(sim[1,:k,:],label = "Global", title = "capital", xlabel = "t", ylabel = "k")

p2 = plot(sim[1,:i,:],label = "Global", title = "investment", xlabel = "t", ylabel = "i")

p3 = plot(sim[1,:d,:],label = "Global", title = "dividends", xlabel = "t", ylabel = "d")

plot(p1,p2,p3,p4, layout = (2,2))   


