using Dolo
using Plots
#  dans le package manager : ] (pour passer en package manager) activate .
model  = yaml_import("examples/models/rbc.yaml")
typeof(model)



sol = perturb(model) # c'est une fonction (une règle de décision) qui dépend de deux arguments
sol = time_iteration(model)
dr = sol.dr
dr([0.0],[3.0])

tab = tabulate(model,dr,:k)

i = tab[V=:i]

k = tab[V=:k]

n = tab[V=:n]

p1 = plot(k,i,xlabel = "k",ylabel="i")
p2 = plot(k,n,xlabel = "k",ylabel="n")

plot(p1,p2,layout=(2,1))
