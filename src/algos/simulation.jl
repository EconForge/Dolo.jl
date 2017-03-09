# function QuantEcon.simulate_2(m::AbstractNumericModel, dr::AbstractDecisionRule,
#                            s0::AbstractVector=m.calibration[:states];
#                            n_exp::Int=0, horizon::Int=40,
#                            seed::Int=42,
#                            forcing_shocks::AbstractMatrix=zeros(0, 0))
#
function simulation(model::AbstractNumericModel,
                           n_exp::Int=0, horizon::Int=40,
                           seed::Int=42,
                           forcing_shocks::AbstractMatrix=zeros(0, 0))


    #forcing_shocks=zeros(0, 0)
    #verbose::Bool=true

    n_exp = max(n_exp, 1)

    # extract data from model
    calib = model.calibration
    params = calib[:parameters]
    #ny = length(calib[:auxiliaries])
    #has_aux = ny > 0

    sigma = (model.calibration.flat[:sig_z])^2

    # calculate initial controls using decision rule
    s0=model.calibration[:states]
    @time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)
    x0 = dr(s0)

    # get number of states and controls
    ns = length(s0)
    nx = length(x0)
    nsx = nx+ns

    # TODO: talk to Pablo and make a decision about this
    s_simul = Array(Float64, n_exp, ns, horizon)
    x_simul = Array(Float64, n_exp, nx, horizon)
    for i in 1:n_exp
      s_simul[i, :, 1] = s0
      x_simul[i, :, 1] = x0
    end
    # NOTE: this will be empty if ny is zero. That's ok. Our call to `cat`  #       below will work either way  y_simul = Array(Float64, n_exp, ny, horizon)


    #using Distributions
    #
    sigma = ones(1,1)*sigma
    n_m = size(sigma,1)
    d = MvNormal(zeros(n_m),sigma)
    epsilons = Array(Float64, n_exp, n_m, horizon)
    for i =1:n_exp
        epsilons[i,:,:] = rand(d,horizon)
    end


    for t in 1:horizon
        #if irf
        #  if !isempty(forcing_shocks) && t < size(forcing_shocks, 2)
        #      epsilons = forcing_shocks[t, :]'
        #  else
        #      epsilons = zeros(size(sigma, 1), 1)
        #  end
        ##else
        ##    epsilons = rand(ϵ_dist, n_exp)'
        #end

        #s = view(s_simul, :, :, t)
        s = copy(view(s_simul, :, :, t))
        #x = view(x_simul, :, :, t)
        x = dr(s)
        # this won't work with s = view(s_simul, :, :, t) even if using vec
        x_simul[:, :, t] = x
        m = view(epsilons,:,:,t)

        if t < horizon
          M = view(epsilons,:,:,t+1)
          ss = view(s_simul, :, :, t+1)
          ss = Dolo.transition!(model, vec(ss), vec(m), vec(s), vec(x), vec(M), params)
          s_simul[:, :, t+1] = ss
          #verbose && @printf "%-8s%-10s%-10s%-10s%-5s\n"  t round(s[1],2) round(s[2],2) round(x[1],2) round(x[2],2)
        end


    end



    cat(2, s_simul, x_simul)::Array{Float64,3}


end
