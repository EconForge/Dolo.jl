import Dolo
import Dolo: Point,ListOfPoints
using StaticArrays


model = Dolo.yaml_import(joinpath(Dolo.pkg_path, "examples","models","rbc_dtcc_iid.yaml"))

grid = model.grid
dprocess = Dolo.discretize(model.exogenous)
n_m = Dolo.n_nodes(dprocess)
init_dr = Dolo.ConstantDecisionRule(model.calibration[:controls])
x0 = [init_dr(i,Dolo.nodes(grid)) for i=1:1]

@time Dolo.improved_time_iteration(model);
@profile Dolo.improved_time_iteration(model);

#
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


R_i = [to_LOP(res[i,:,:]) for i=1:size(res,1)]
D_i =  [to_LOJ(dres[i,:,:,:]) for i=1:size(res,1)]
J_ij = typeof(D_i[1])[to_LOJ(jres[i,j,:,:,:]) for i=1:size(jres,1), j=1:size(jres,2)]
S_ij = typeof(R_i[1])[to_LOP(fut_S[i,j,:,:]) for i=1:size(fut_S,1), j=1:size(fut_S,2)]

ddr_filt = Dolo.CachedDecisionRule(dprocess, grid, R_i)

Π_i = deepcopy(R_i)
π_i = deepcopy(R_i)
Dinv = deepcopy(D_i)
for i=1:length(Dinv)
    invert!(Dinv[i])
end
M_ij = deepcopy(J_ij)
for i=1:size(M_ij,1)
    for j=1:size(M_ij,2)
        premult!(Dinv[i],M_ij[i,j])
        M_ij[i,j][:] *= Dolo.iweight(dprocess,i,j)
    end
end



module NewModule

    import Dolo
    import Dolo: Point,ListOfPoints
    using StaticArrays

    function to_LOP(mat::Matrix{Float64})
        N,d = size(mat)
        reinterpret(Point{d}, mat', (N,))
    end

    function to_LOJ(mat::Array{Float64})
        # list of jacobians
        N,d = size(mat)
        reinterpret(SMatrix{d,d,Float64,d*d}, permutedims(mat,[2,3,1]), (N,))
    end

    function invert!(A)
        # A[i] <- (A[i])^(-1)
        N = length(A)
        for n=1:N
            A[n] = inv(A[n])
        end
    end

    function premult!(A,B)
        # B[i] <- A[i]*B[i]
        N = length(A)
        for n=1:N
            B[n] = A[n]*B[n]
        end
    end

    function addmul!(O,A,B)
        # O[i] <- A[i]*B[i]
        N = length(A)
        for n=1:N
            O[n] += A[n]*B[n]
        end
    end


    import Base


    function Base.maxabs(v::ListOfPoints{d}) where d
        m = 0.0
        for n = 1:length(v)
            mm = maximum(abs, v[n])
            if mm>m
                m = mm
            end
        end
        m
    end

    function Base.maxabs(v::Vector{ListOfPoints{d}}) where d
        maximum(maxabs.(v))
    end




    function d_filt_dx!(Π_i::Vector{ListOfPoints{n_x}}, π_i::Vector{ListOfPoints{n_x}}, M_ij, S_ij::Matrix{ListOfPoints{d}}, dumdr) where n_x where d

      n_m,n_mt = size(M_ij)
      Dolo.set_values!(dumdr,π_i)
      for i in 1:n_m
        Π_i[i][:] *= 0
        for j in 1:n_mt
          A = M_ij[i,j]
          B = dumdr(i,j,S_ij[i,j])
          addmul!(Π_i[i],A,B)
        end
      end

    end

    function do_it_many_times(R_i, M_ij, S_ij, ddr_filt; N=10)
        Π_i = deepcopy(R_i)
        π_i = deepcopy(R_i)
        errors = [maxabs(R_i)]
        for it=1:N
            d_filt_dx!(Π_i, π_i, M_ij, S_ij, ddr_filt)
            π_i = Π_i
            push!(errors, maxabs(Π_i))
        end
        return errors, π_i
    end

end


import NewModule

@time errs, sol = do_it_many_times(R_i, M_ij, S_ij, ddr_filt; N=2)
@time errs, sol = NewModule.do_it_many_times(R_i, M_ij, S_ij, ddr_filt; N=476)


rats = errs[2:end]./errs[1:(end-1)]

size(R_i[1])
