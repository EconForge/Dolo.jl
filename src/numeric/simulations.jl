function QuantEcon.simulate(m::AbstractNumericModel, dr::AbstractDecisionRule,
                            s0::AbstractVector=m.calibration[:states];
                            n_exp::Int=0, horizon::Int=40,
                            seed::Int=42,
                            forcing_shocks::AbstractMatrix=zeros(0, 0))

    irf = n_exp == 0 ? true : false
    n_exp = max(n_exp, 1)

    # extract data from model
    calib = m.calibration
    params = calib[:parameters]
    ny = length(calib[:auxiliaries])
    has_aux = ny > 0

    sigma = m.options.distribution.sigma

    # calculate initial controls using decision rule
    x0 = dr(s0)

    # get number of states and controls
    ns = length(s0)
    nx = length(x0)

    # TODO: talk to Pablo and make a decision about this
    s_simul = Array(Float64, n_exp, ns, horizon)
    x_simul = Array(Float64, n_exp, nx, horizon)
    for i in 1:n_exp
        s_simul[i, :, 1] = s0
        x_simul[i, :, 1] = x0
    end

    # NOTE: this will be empty if ny is zero. That's ok. Our call to `cat`
    #       below will work either way
    y_simul = Array(Float64, n_exp, ny, horizon)

    # extract functions that we'll use
    f = m.functions.arbitrage
    g = m.functions.transition
    a = m.functions.auxiliary

    srand(seed)
    ϵ_dist = MvNormal(zeros(size(sigma, 1)), sigma)

    for t in 1:horizon
        if irf
            if !isempty(forcing_shocks) && t < size(forcing_shocks, 2)
                epsilons = forcing_shocks[t, :]'
            else
                epsilons = zeros(size(sigma, 1), 1)
            end
        else
            epsilons = rand(ϵ_dist, n_exp)'
        end

        s = view(s_simul, :, :, t)

        # extract x and then fill it
        x = view(x_simul, :, :, t)
        dr(s, x)

        if has_aux
            # extract a and then fill it
            y = view(y_simul, :, :, t)
            evaluate!(a, s, x, params, y)
        end

        if t < horizon
            ss = view(s_simul, :, :, t+1)
            evaluate!(g, s, x, epsilons, params, ss)
        end
    end

    # now stitch everything together
    cat(2, s_simul, x_simul, y_simul)::Array{Float64,3}

end


#=
url = "https://raw.githubusercontent.com/EconForge/dolo/c8bd2e3f2f5402f687beb7949d49deefda6a5fc6/examples/models/rbc.yaml"
m = yaml_import(url)
dr = linear_solve(m)
sim = simulate(m, dr)
s0 = m.calibration[:states]
horizon = 40
n_exp = 4
seed = 42
discard=false
forcing_shocks = zeros(0, 0)

@code_warntype simulate(m, dr, s0, 0, 40, 42, zeros(0, 0))
=#
