using Dolo

model = yaml_import("examples/models/rbc.yaml")


# sol = Dolo.improved_time_iteration(model)
@time res = Dolo.time_iteration(model, ignore_constraints=true, verbose=false);


F = Dolo.Euler(model)


J = Dolo.df_A(F, F.x0, F.x0)

Jii = Dolo.df_A(F, F.x0, F.x0, inplace=true)

z0 = F.x0


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


J = Dolo.DiffFun!(fun!, F.x0, 1e-8, out)