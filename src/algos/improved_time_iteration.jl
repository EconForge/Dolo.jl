include("ITI_additional.jl")
using JLD

"""
Computes a global solution for a model via backward Improved Time Iteration. The algorithm is applied to the residuals of the arbitrage equations. The idea is to solve the system G(x) = 0 as a big nonlinear system in x, where the inverted Jacobian matrix is approximated by an infinite sum (Neumann series).

If the initial guess for the decision rule is not explicitly provided, the initial guess is provided by `ConstantDecisionRule`.
If the stochastic process for the model is not explicitly provided, the process is taken from the default provided by the model object, `model.exogenous`

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `dprocess`: The stochastic process associated with the exogenous variables in the model.
* `init_dr`: Initial guess for the decision rule.
* `maxbsteps` Maximum number of backsteps.
* `verbose` Set "true" if you would like to see the details of the infinite sum convergence.
* `smaxit` Maximum number of iterations to compute the Neumann series.
* `complementarities`
* `compute_radius`
* `details` If false returns only a decision rule dr
# Returns
* `dr`: Solved decision rule.
* `details` about the iterations is specified.
"""

function improved_time_iteration(model::AbstractModel, dprocess::AbstractDiscretizedProcess,
                                 init_dr::AbstractDecisionRule, grid;
                                 maxbsteps::Int=10, verbose::Bool=true, verbose_jac::Bool=false,
                                 tol::Float64=1e-8, smaxit::Int=500, maxit::Int=1000,
                                 complementarities::Bool=false, compute_radius::Bool=false, details::Bool=true)

   # x_lb = model.functions['controls_lb']
   # x_ub = model.functions['controls_ub']

   parms = model.calibration[:parameters]

  #  n_m = Dolo.n_nodes(dprocess) # number of exo states today
   n_m = max(n_nodes(dprocess), 1) # number of exo states today
   n_mt = n_inodes(dprocess,1)  # number of exo states tomorrow
   n_s = length(model.symbols[:states]) # number of endo states

   s = nodes(grid)
   N_s = size(s,1)
   n_x = size(model.calibration[:controls],1)
  #  N_m = Dolo.n_nodes(dprocess) # number of grid points for exo_vars

  #  x0 = [repmat(model.calibration[:controls]',N_s) for i in 1:N_m] #n_x N_s n_m
   x0 = [init_dr(i, s) for i=1:n_m]
   ddr=CachedDecisionRule(dprocess, grid, x0)
   ddr_filt = CachedDecisionRule(dprocess, grid, x0)
   set_values!(ddr,x0)

   steps = 0.5.^collect(0:maxbsteps)

   if complementarities == true
     x_lb = Array{Float64,2}[cat(1, [controls_lb(model, node(dprocess, i), s[n, :], parms)' for n=1:N_s]...) for i=1:n_m]
     x_ub = Array{Float64,2}[cat(1, [controls_ub(model, node(dprocess, i), s[n, :], parms)' for n=1:N_s]...) for i=1:n_m]
   end


   x=x0
   ## memory allocation
   jres = zeros(n_m,n_mt,N_s,n_x,n_x)
   S_ij = zeros(n_m,n_mt,N_s,n_s)

   ######### Loop     for it in range(maxit):
   it=0
   it_invert=0
   res_init = euler_residuals(model,s,x,ddr,dprocess,parms ,set_dr=false, jres=jres, S_ij=S_ij)
   err_0 = abs(maximum(res_init))
   err_2= err_0
   lam0=0.0

   verbose && println(repeat("-", 120))
   verbose && println("N\tf_x\t\td_x\tTime_residuals\tTime_inversion\tTime_search\tLambda_0\tN_invert\tN_search\t")
   verbose && println(repeat("-", 120))


   if compute_radius == true
     res=zeros(res_init)
     dres = zeros(N_s*n_m, n_x, n_x)
   end

   while it <= maxit && err_0>tol
      it += 1

      jres = zeros(n_m,n_mt,N_s,n_x,n_x)
      S_ij = zeros(n_m,n_mt,N_s,n_s)

      t1 = time();

      # compute derivatives and residuals:
      # res: residuals
      # dres: derivatives w.r.t. x
      # jres: derivatives w.r.t. ~x
      # fut_S: future states
      set_values!(ddr,x)

      ff = SerialDifferentiableFunction(u-> euler_residuals(model, s, u,ddr,dprocess,parms;
                                        with_jres=false,set_dr=false))

      res, dres = ff(x)

      # dres = permutedims(dres, [axisdim(dres, Axis{:n_v}),axisdim(dres, Axis{:N}),axisdim(dres, Axis{:n_x})])
      dres = reshape(dres, n_m, N_s, n_x, n_x)
      junk, jres, fut_S = euler_residuals(model, s, x,ddr,dprocess,parms, with_jres=true,set_dr=false, jres=jres, S_ij=S_ij)

      if complementarities == true
        for i_ms in 1:n_m
           dx =  x[i_ms] - x_lb[i_ms]
          #  res.data[i_ms,:,:], dres[i_ms,:,:,:], jres[i_ms,:,:,:,:] = smooth_right(res.data[i_ms,:,:], dres[i_ms,:,:,:], jres[i_ms,:,:,:,:], dx)
           res[i_ms,:,:], dres[i_ms,:,:,:], jres[i_ms,:,:,:,:] = smooth_right(res[i_ms,:,:], dres[i_ms,:,:,:], jres[i_ms,:,:,:,:], dx)
        end

        res *= -1
        dres *= -1
        jres *= -1

        # i_ms=1
        # smooth_right(res.data[i_ms,:,:], dres[i_ms,:,:,:], jres[i_ms,:,:,:,:], dx; pos = -1.0)
        for i_ms in 1:n_m
           dx =  x_ub[i_ms] -x[i_ms]
           res[i_ms,:,:], dres[i_ms,:,:,:], jres[i_ms,:,:,:,:] = smooth_right(res[i_ms,:,:], dres[i_ms,:,:,:], jres[i_ms,:,:,:,:], dx; pos = -1.0)
        end
      end

      err_0 = abs(maximum(res))

      jres *= -1.0
      M=jres

      # X=zeros(n_m,N_s,n_x,n_x)
      for i_m in 1:n_m
        X = copy(dres[i_m,:,:,:])
        for j_m in 1:n_mt
            for n in 1:N_s
                M[i_m,j_m,n,:,:] = X[n,:,:]\M[i_m,j_m,n,:,:]
            end
        end
      end

    #   if it==1
    #       save("myfile.jld", "res", res, "dres", dres, "jres", jres, "fut_S", fut_S)
    #   end

      ####################
      # Invert Jacobians
      t2 = time();
      tot, it_invert, lam0, errors = invert_jac(res,dres,jres,fut_S, ddr_filt; verbose=verbose_jac, maxit = smaxit)

      t3 = time();

      i_bckstps=0
      new_err=err_0
      new_x = x
      while new_err>=err_0 && i_bckstps<length(steps)
        i_bckstps +=1
        new_x = x-destack0(tot, n_m)*steps[i_bckstps]
        new_res = euler_residuals(model, s, new_x,ddr,dprocess,parms,set_dr=true)

        if complementarities == true
          for i_ms in 1:n_m
             dx =  new_x[i_ms]-x_lb[i_ms]
             new_res[i_ms,:,:] = smooth_right(new_res[i_ms,:,:], dx)
          end
          for i_ms in 1:n_m
             dx =  x_ub[i_ms] - new_x[i_ms]
             new_res[i_ms,:,:] = smooth_right(-new_res[i_ms,:,:], dx)
          end
        end

        new_err = maximum(abs, new_res)
      end
      err_2 = maximum(abs,tot)

      t4 = time();

      x = new_x
      verbose && @printf "%-6i% -10e% -17e% -15.4f% -15.4f% -15.5f% -17.3f%-17i%-5i\n" it  err_0  err_2  t2-t1 t3-t2 t4-t3 lam0 it_invert i_bckstps

   end
   verbose && println(repeat("-", 120))
   set_values!(ddr,x)

   if compute_radius == true
     lam, lam_max, lambdas = radius_jac(res,dres,jres,S_ij,ddr_filt)
   end

   if !details
     return ddr.dr
   else
     converged = err_0<tol
     if !compute_radius
       return ImprovedTimeIterationResult(ddr.dr, it, err_0, err_2, converged, complementarities, tol, lam0, it_invert, 5.0)
     else
       return ImprovedTimeIterationResult(ddr.dr, it, err_0, err_2, converged, tol, lam0, it_invert, 5.0), (lam, lam_max, lambdas)
     end
   end

end

function improved_time_iteration(model:: AbstractModel, dprocess::AbstractDiscretizedProcess,
                                 init_dr::AbstractDecisionRule;grid=Dict(), kwargs...)
    grid = get_grid(model, options=grid)
    return improved_time_iteration(model, dprocess, init_dr, grid;  kwargs...)
end

function improved_time_iteration(model, dprocess::AbstractDiscretizedProcess; grid=Dict(), kwargs...)

    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return improved_time_iteration(model, dprocess, init_dr; grid=grid, kwargs...)
end


function  improved_time_iteration(model, init_dr; grid=Dict(), kwargs...)
    dprocess = discretize( model.exogenous )
    return improved_time_iteration(model, dprocess, init_dr; grid=grid, kwargs...)
end

# function improved_time_iteration(model, maxbsteps::Int=10, verbose::Bool=false,
#                                  tol::Float64=1e-8, smaxit::Int=500, maxit::Int=1000,
#                                  complementarities::Bool=true, compute_radius::Bool=false)
#     dprocess = Dolo.discretize( model.exogenous )
#     init_dr = Dolo.ConstantDecisionRule(model.calibration[:controls])
#     return improved_time_iteration(model, dprocess, init_dr, maxbsteps, verbose,tol,
#                                    smaxit, maxit,complementarities, compute_radius)
# end

function improved_time_iteration(model; grid=Dict(), kwargs...)
    dprocess = discretize( model.exogenous )
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return improved_time_iteration(model, dprocess, init_dr; grid=grid, kwargs...)
end
