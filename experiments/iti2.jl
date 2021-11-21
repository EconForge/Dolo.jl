using Dolo

model = yaml_import("examples/models/rbc.yaml")


sol = Dolo.improved_time_iteration(model; ignore_constraints=false)



@time res = Dolo.time_iteration(model, ignore_constraints=false, verbose=false);


F = Dolo.Euler(model)


@time F(F.x0, F.x0);

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