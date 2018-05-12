import Dolo
import Dolo: Point,ListOfPoints
using StaticArrays


model = Dolo.yaml_import(joinpath(Dolo.pkg_path, "examples","models","rbc.yaml"));


@time sol_iti = Dolo.improved_time_iteration(model; verbose=true);
@time sol_iti = Dolo.improved_time_iteration(model; verbose=true);


@time sol_ti = Dolo.time_iteration(model, verbose=true; tol_η=1e-10, maxit=2000)
exit()
# @profile Dolo.improved_time_iteration(model);
#
# #
# f = open("out.jll", "r")
# out = deserialize(f)
# close(f)
#
#
# res, dres, jres, fut_S = out

using JLD
ss = load("myfile.jld")
res = ss["res"]
dres = ss["dres"]
jres = ss["jres"]
fut_S = ss["fut_S"]



module NewModule

    import Dolo
    import Dolo: Point,ListOfPoints
    using StaticArrays


end



import NewModule

R_i, D_i, J_ij, S_ij, M_ij = NewModule.reorder_data(dprocess, res, dres, jres, fut_S)

ddr_filt = Dolo.CachedDecisionRule(dprocess, grid, R_i)


@time tot0, n0, lam0 = Dolo.invert_jac(res,dres,jres,fut_S,ddr_filt)

n0

@time tot, n, lam, errors= NewModule.invert_jac(res,dres,jres,fut_S,ddr_filt)



print("Invert me many times")

function invert_me_1(N, res,dres,jres,fut_S,ddr_filt;maxit=1000)
    tt = 0
    for n=1:N
        tot0, n0, lam0 = Dolo.invert_jac(res,dres,jres,fut_S,ddr_filt;maxit=maxit)
        tt += n0
    end
    return tt
end

function invert_me_2(N, res,dres,jres,fut_S,ddr_filt;maxit=1000)
    tt = 0
    for n=1:N
        tot0, n0, lam0 = NewModule.invert_jac(res,dres,jres,fut_S,ddr_filt;maxit=maxit)
        tt += n0
    end
    return tt
end

@time invert_me_1(1,res,dres,jres,fut_S,ddr_filt;maxit=2)
@time invert_me_2(1,res,dres,jres,fut_S,ddr_filt;maxit=2)




@time invert_me_1(1,res,dres,jres,fut_S,ddr_filt;maxit=1000)
@time invert_me_2(1,res,dres,jres,fut_S,ddr_filt;maxit=1000)
