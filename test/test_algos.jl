path = Pkg.dir("Dolo")

import Dolo
import YAML

@testset "testing model algos" begin

    @testset "testing mc models" begin
        fn = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
        # fn = joinpath(path,"examples","models","LAMP.yaml")
        model_mc = Dolo.yaml_import(fn)

        drc = Dolo.ConstantDecisionRule(model_mc.calibration[:controls])
        # @time dr0, drv0 = Dolo.solve_policy(model_mc, drc) #, verbose=true, maxit=10000 )
        @time dr = Dolo.time_iteration(model_mc, verbose=false, maxit=10000)
        @time drv = Dolo.evaluate_policy(model_mc, dr, verbose=false) #, verbose=true, maxit=10000)
        @time drd = Dolo.time_iteration_direct(model_mc, dr, verbose=false) #, maxit=500)

        sim = Dolo.simulate(model_mc, dr) #; N=100, T=20)

    end

# # compare with prerecorded values
# kvec = linspace(dr.grid.min[1],dr.grid.max[1],10)
# nvec = [dr(1,[k])[1] for k in kvec]
# ivec = [dr(1,[k])[2] for k in kvec]
# # compare  time_iteration_direct
# nvec_d = [drd(1,[k])[1] for k in kvec]
# ivec_d = [drd(1,[k])[2] for k in kvec]
# @assert maximum(abs, nvec_d-nvec)<1e-4
#
# # compare  vfi
# nvec_0 = [dr0(1,[k])[1] for k in kvec]
# ivec_0 = [dr0(1,[k])[2] for k in kvec]
# @assert maximum(abs, nvec_0-nvec)<1e-4

# let's redo when model is stable !
# ivec_test = [0.295977,  0.257538,  0.21566,  0.173564,  0.132103,  0.0915598,  0.0520067,  0.0134661,  7.01983e-6, 3.40994e-17]
# nvec_test = [ 0.391997,  0.348033,  0.318369,  0.296276,  0.278821,  0.264487,  0.25239 ,  0.241974,  0.236604,  0.233779 ]
# @assert maximum(abs(ivec-ivec_test))<1e-5
# @assert maximum(abs(nvec-nvec_test))<1e-5

    @testset "testing iid models" begin


        import Dolo
        path = Pkg.dir("Dolo")
        fn = joinpath(path,"examples","models","rbc_dtcc_iid.yaml")
        model = Dolo.yaml_import(fn)

        @time dr = Dolo.perturbate(model)

        drc = Dolo.ConstantDecisionRule(model.calibration[:controls])

        @time dr = Dolo.time_iteration(model, maxit=100, verbose=false)
        @time res = Dolo.time_iteration(model, dr; maxit=100, details=true)

        @time dr0, drv0 = Dolo.solve_policy(model, drc) #;, verbose=true, maxit=1000 )
        @time res = Dolo.solve_policy(model, drc; details=true) #;, verbose=true, maxit=1000 )

        @time drd = Dolo.time_iteration_direct(model) #, maxit=1000, verbose=true)

        @time dr = Dolo.time_iteration_direct(model, drd) #, maxit=500, verbose=true)
        @time res = Dolo.time_iteration_direct(model, drc; details=true)
        @time drv = Dolo.evaluate_policy(model, dr; verbose=false)

        #
        Dolo.simulate(model, dr)

        s0 = model.calibration[:states]+0.1
        sim = Dolo.simulate(model, dr, s0)
        Dolo.simulate(model, dr; N=10)


        res = Dolo.response(model.exogenous, [0.01])


        irf = Dolo.response(model, dr, :e_z)

        irf[:k]
    end

# kvec = linspace(dr.grid.min[1],dr.grid.max[1],10)
# nvec = [dr(1,[k])[1] for k in kvec]
# ivec = [dr(1,[k])[2] for k in kvec]
# nvec_d = [drd(1,[k])[1] for k in kvec]
# ivec_d = [drd(1,[k])[2] for k in kvec]
# nvec_0 = [dr0(1,[k])[1] for k in kvec]
# ivec_0 = [dr0(1,[k])[2] for k in kvec]
#
# @assert maximum(abs, nvec_d-nvec)<1e-5
# @assert maximum(abs, nvec_0-nvec)<1e-5 # not satisfied right now (see tol. of optimizer)

    @testset "testing ar1 models" begin

        # AR1 model: this one should be exactly equivalent to rbc_dtcc_ar1
        import Dolo
        fn = Pkg.dir("Dolo","examples","models","rbc_dtcc_ar1.yaml")
        model = Dolo.yaml_import(fn)
        dp = Dolo.discretize(model.exogenous)

        @time dr = Dolo.perturbate(model)
        cdr =Dolo.CachedDecisionRule(dr,dp)
        @time dr = Dolo.time_iteration(model, model.exogenous, cdr)
        @time dr = Dolo.time_iteration(model)
        @time dr0, drv0 = Dolo.solve_policy(model, cdr, verbose=false) #, maxit=10000 )
        @time drv = Dolo.evaluate_policy(model, dr, verbose=false, maxit=10000)
        @time drd = Dolo.time_iteration_direct(model, dr) #, verbose=true) #, maxit=500)
        #
        model.symbols[:exogenous]

        Dolo.simulate(model, dr, N=10)

        Dolo.response(model.exogenous, [0.01])

        sim = Dolo.response(model, dr, :z)

        sim[:z]
    end

end
