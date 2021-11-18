using Dolo


# model = Model("examples/models/rbc.yaml")
model = Model("examples/models/consumption_savings_iid.yaml")

F = Dolo.Euler(model, ignore_constraints=false);


F(F.x0, F.x0)

r = F(F.x0, F.x0; set_future=false)
j = Dolo.df_A(F, F.x0, F.x0)

import Dolo: norm

import Dolo: MSM
import Dolo: df_A

f = u->F(u, F.x0; set_future=false)
df = u->df_A(F, u, F.x0; set_future=false)

function newton(fun, dfun, x0::MSM; maxit=50, dampen=1.0, verbose=false)

    r0 = fun(x0)

    err = norm(r0)

    local x1

    bcksteps = 5

    for n=1:maxit

        r0 = fun(x0)
        err_0 = norm(r0)

        j = dfun(x0)
        δ = j\r0

        for i=0:(bcksteps-1)
            u = 0.5^i
            x1 = x0 - δ*u
            r1 = fun(x1)
            err_1 = norm(r1)
            if verbose
                println( "-    $i: ", err_1)
            end
            if err_1<err_0
                break
            end
        end

        if verbose
            println(n, " | ", norm(r0), " | ",  norm(δ))
        end

        x0 = x1

    end

    return x0

end

sol = newton(f, df, F.x0; maxit=20, verbose=true)

norm(F(sol, F.x0))


F(sol, sol; set_future=true);


function loop(F, K)

    x0 = F.x0
    x1 = x0

    for k=1:K

        res = F(x0, x0; set_future=true)
        ε = norm(res)

        f = u->F(u, x1; set_future=false)
        df = u->df_A(F, u, x1; set_future=false)

        x1 = newton(f, df, x0; verbose=false)

        δ = x1 - x0
        η = norm(δ)

        println(ε, " : ", η)

        x0 = x1
    end

    return x1

end


loop(F, 10);


# tab = Dolo.tabulate(model, F.dr.dr, :k)
# plot(tab[:k], tab[:k])
# plot!(tab[:j], tab[:c])



using SimplePlots


tab = Dolo.tabulate(model, F.dr.dr, :w)
plot(tab[:w], tab[:w])
plot!(tab[:w], tab[:c])


# l = Dolo.df_B(F, F.x0, F.x0)



# Dolo.newton(u->F(u,F.x0), F.x0)



sol = Dolo.time_iteration(model, maxit=1000);
tab = tabulate(model, sol.dr, :m)
pl = plot(tab[:m], tab[:m])
plot!(pl, tab[:m], tab[:c])



Dolo.time_iteration(model; ignore_constraints=true);


Dolo.improved_time_iteration(model);


# grid, exo = Dolo.discretize(model)

println("Endo Grid size: $(Dolo.n_nodes(grid.endo))")
println("Exo Grid size: $(Dolo.n_nodes(grid.exo))")
println("Exo Grid I size: $(Dolo.n_inodes(exo, 1))")


grid, exo = Dolo.discretize(model;)

println("Endo Grid size: $(Dolo.n_nodes(grid.endo))")
println("Exo Grid size: $(Dolo.n_nodes(grid.exo))")
println("Exo Grid I size: $(Dolo.n_inodes(exo, 1))")

sol0 = time_iteration(model)



F = Dolo.Euler(model)

@time Dolo.loop_iti(F, F.x0)

@time Dolo.loop_ti(F, F.x0)


# μ = Dolo.ergodic_distribution(model, sol)



# using Plots

