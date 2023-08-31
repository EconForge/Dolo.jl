

using Dolo
const NL = Dolo

root_dir = pkgdir(Dolo)
# model = yaml_import("$(root_dir)/examples/ymodels/rbc_mc.yaml")
model = include("$(root_dir)/examples/ymodels/rbc_mc.jl")

model.calibration

model2 = NL.recalibrate(model; β=0.8)

dmodel = NL.discretize(model)

sol = NL.time_iteration(dmodel)

tab = tabulate(model, sol.dr, :k;)

sim = NL.simulate(model, sol.dr)



model = include("$(root_dir)/examples/ymodels/rbc_iid.jl")

s0 = NL.draw(model.states)

x = SVector( 0.2, 0.3)

gen = NL.τ(model, s0, x)

NL.calibrated(NL.QP, model, :states)
# NL.calibrated(NL.QP, model)

NL.QP(model.states, [0.2, 0.4])

    

NL.QP(model.states, SVector(0.0, 5.5))


@time sol = Dolo.time_iteration(model);
@time Dolo.time_iteration(model, improve=true, improve_wait=0, improve_K=500);


# here the problem is the DFun is not initialized with the right
# variables

tab = tabulate(model, sol.dr, :k)

s0 = NL.draw(model.states)

sim = NL.simulate(model, sol.dr, s0)



# model = yaml_import("$(root_dir)/examples/ymodels/consumption_savings_iid.yaml")

# # s0 = [NL.enum(dm.grid)...][50]
# # x = SVector( 0.5)
# # function f(dm, s0, x)
# #     r = sum( w*S.val for (w,S) in NL.τ(dm, s0, x))
# #     return sum(r)
# # end
# # @time f(dm, s0, x)
# # @code_warntype f(dm, s0, x)


# function initial_guess(model::typeof(model), s::NamedTuple)
#     c = s.w*0.9
#     return (;c)
# end

# @time Dolo.time_iteration(model; verbose=true, trace=true);
# @time Dolo.time_iteration(model; improve=true, improve_K=500);

# @time Dolo.vfi(model; improve=true);



# model = yaml_import("$(root_dir)/examples/ymodels/consumption_savings.yaml")

# # dm = NL.discretize(model)


# # function initial_guess(model::typeof(model), s::NamedTuple)
# #     c = s.w*0.9
# #     return (;c)
# # end

# # s0 = [NL.enum(dm.grid)...][50]

# # x = SVector( 0.5)

# # function f(dm, s0, x)
# #     r = sum( w*S.val for (w,S) in NL.τ(dm, s0, x))
# #     return sum(r)
# # end

# # @time f(dm, s0, x)

# # @code_warntype f(dm, s0, x)



# # dm = Dolo.discretize(model)

# @time Dolo.time_iteration(model; verbose=true, T=10);
# @time Dolo.time_iteration(model; improve=true, improve_K=1000);

# @time Dolo.vfi(model; improve=true);



model2 = include("$(root_dir)/examples/ymodels/consumption_savings.jl")

# # problem: 
# dm2 = Dolo.discretize(model2, Dict(:endo=>Dict(:n=>100), :exo=>Dict(:n=>7)))

# s1 = [NL.enum(dm2.grid)...][50]
# @time f(dm2, s1, x)
# @code_warntype f(dm2, s1, x)
# @time sum( w*S.val for (w,S) in NL.τ(dm2, s1, x))


# @time Dolo.time_iteration(dm2, verbose=true, T=10);

@time Dolo.time_iteration(model; verbose=true);
@time Dolo.time_iteration(model; improve=true, improve_K=1000);

# @time Dolo.vfi(model; improve=true);





model = include("$(root_dir)/examples/ymodels/rbc_ar1.jl")

# # problem: 
# dm = Dolo.discretize(model)
# x = SVector(0.1, 0.4)

# s0 = [NL.enum(dm.grid)...][50]
# @time f(dm, s0, x)
# @code_warntype f(dm, s0, x)
# @time sum( w*S.val for (w,S) in NL.τ(dm, s0, x))


# @time Dolo.time_iteration(dm2; verbose=true, T=10);

@time Dolo.time_iteration(model; verbose=true);
@time Dolo.time_iteration(model; improve=true, improve_K=1000);

# @time Dolo.vfi(model; improve=true);
