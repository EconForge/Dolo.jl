using Dolo

model = yaml_import("examples/models/consumption_savings_iid.yaml")


model = yaml_import("examples/models/rbc_mc.yaml")



sol = Dolo.time_iteration(model, verbose=true, ignore_constraints=false, maxit=5);


simulate(model, sol.dr)


@time sol = Dolo.improved_time_iteration(model; verbose=false, ignore_constraints=false);


using SimplePlots

tab = tabulate(model, sol.dr, :w)

plot(tab[:w], tab[:c])

function step(model)
    sol = Dolo.improved_time_iteration(model; verbose=false, ignore_constraints=false);
    dist = Dolo.ergodic_distribution(model, sol);
end


@time step(model);

using Dolo: Euler, MSM
using StaticArrays

function dFdp(F::Euler, x0::MSM, x1::MSM, p0::SVector, dp::SVector, ε=1e-8)
    
    no = maximum(abs, dp)
    println(no)
    p1 = p0 + (dp/no)*ε

    f0 = F(x0, x1; p=p0)
    f1 = F(x0, x1; p=p1)

    df = ((f1-f0)/ε) * no

    return df

end


function update_guess(F, x0, p0, dp, tol=1e-8, maxit=1000)
    df = dFdp(F, x0, x0, p0, dp)
    J = Dolo.df_A(F, x0, x0)
    L = Dolo.df_B(F, x0, x0)
    Dolo.prediv!(L, J) 
    Dolo.mult!(L, -1.0)
    dπ = -J\df
    dx = dπ
    i = 0
    while i<=maxit
        i += 1
        err = Dolo.norm(dπ)
        if err<tol
            break
        end
        dπ = L*dπ
        dx += dπ
    end
    return dx, i

end


sol =  improved_time_iteration(model, verbose=false);


dp = SVector(0.01, 0, 0, 0, 0, 0)






@time res = Dolo.time_iteration(model, ignore_constraints=false, verbose=false);

xx = MSM([sol.dr(1, F.grid_endo.nodes)])

@time dx, i = update_guess(F, xx, F.p, dp);

F(xx+dx, xx+dx; p=p+dp, set_future=true);

@time J = Dolo.df_A(F, F.x0, F.x0, inplace=false);
@time J = Dolo.df_A(F, F.x0, F.x0, inplace=true);

Jii = Dolo.df_A(F, F.x0, F.x0, inplace=true);

f0 = z0*0
fi = z0*0
xi = z0*0

using StaticArrays
using Dolo: MSM
n_x = length(f0.data[1])
N = length(xi.data)
cata = zeros(SMatrix{n_x, n_x, Float64, n_x*n_x}, N)
J = MSM(cata, f0.sizes)
out = (f0, fi, xi, J)


function fun!(out, u)
    r = F(u,F.x0)
    out.data[:] .= r.data
end

fun!(f0, F.x0)

fp = f0*0


fp, Jp = Dolo.DiffFun!(fun!, F.x0, 1e-8, out);

fp.data[1]

F(F.x0, F.x0).data[1]
