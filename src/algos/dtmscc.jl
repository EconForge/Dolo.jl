using QuantEcon: gridmake

using Optim

using Interpolations
using splines # temporary

#####
# helper functions/types (to be moved/removed)
#####

mlinspace(a,b,dims) = gridmake([linspace(a[i],b[i],dims[i]) for i=1:length(dims)]...)


type MarkovChain
    P::Array{Float64, 2}
    Q::Array{Float64, 2}
end

type ApproximationSpace
  a:: Array{Float64, 1}
  b:: Array{Float64, 1}
  dims:: Array{Int64, 1}
end

###
# Object to represent a decision rule/policy which depends on a discrete state
###

type MixedDecisionRule
    n_mc:: Int64
    a:: Array{Float64, 1}
    b:: Array{Float64, 1}
    dims:: Array{Int64, 1}
    values:: Array{Float64, 3}
    interpolants:: Any
    function MixedDecisionRule(n_mc, a, b, dims)
        d = length(dims)
        values = zeros(0,0,0)
        interpolant = Union{}
        return new(n_mc, a, b, dims, values, interpolant)
    end
end

function set_values(mdr::MixedDecisionRule, values)
    n_mc = mdr.n_mc
    dims = mdr.dims
    n_x = size(values)[end]
    d = length(mdr.dims)
    knots = [linspace(mdr.a[i], mdr.b[i], mdr.dims[i]) for i=1:d]
    mdr.interpolants = [ scale(interpolate(copy(reshape(slice(values,i_mc,:,i_x),dims...)), BSpline(Linear()), OnGrid()), knots...) for i_mc in 1:n_mc, i_x in 1:n_x  ]
end

function evaluate(mdr::MixedDecisionRule, i::Int64, s::Array{Float64,1})
    n_x = size(mdr.interpolants,2)
    return Float64[mdr.interpolants[i,j][s...] for j=1:n_x ]
end

function evaluate(mdr::MixedDecisionRule, i::Int64, s::Array{Float64,2})
    return cat(1, [evaluate(mdr, i, copy(slice(s,n,:)))' for n=1:size(s,1)]... )
end

####
# residuals
####

function residuals(model::DTMSCCModel, calibration::ModelCalibration)
    m = calibration[:markov_states]
    s = calibration[:states]
    x = calibration[:controls]
    p = calibration[:parameters]
    res = Dict()
    res[:transition] = (s-evaluate(model.functions.transition,m,s,x,m,p))
    res[:arbitrage] = (evaluate(model.functions.arbitrage,m,s,x,m,s,x,p))
    return res
end

residuals(model::DTMSCCModel) = residuals(model, model.calibration)

####
# construction of initial guess
####

type MixedConstantPolicy
    x:: Array{Float64,1}
end

evaluate(mdr::MixedConstantPolicy, i, s::AbstractArray{Float64,1}) = mdr.x
evaluate(mdr::MixedConstantPolicy, i, s::AbstractArray{Float64,2}) = repmat(mdr.x',size(s,1),1)
#
constant_guess(model::DTMSCCModel) = MixedConstantPolicy(model.calibration[:controls])

###
# simple time iteration
###

function step_residuals(f, g, s::Array{Float64,2}, controls::Array{Float64,3}, p::Array{Float64,1}, mdr, P::Array{Float64,2}, Q::Array{Float64,2})

    N = (size(s, 1))
    n_mc = size(P,1)
    # N == size(x,2)
    n_x = size(controls,3)
    # n_m = size(P,2)
    # size(P,1)==size(P,2)==size(Q,1)

    res = zeros(size(controls))

    # for each discrete state today
    for i_mc=1:n_mc
        # *vector* of markov states today
        m = slice(P,i_mc,:)
        # *matrix* of controls taken today
        x = slice(controls,i_mc,:,:)
        # for each discrete state tomorrow
        for j_mc=1:n_mc
            prob = Q[i_mc, j_mc]
            # *vector* of markov states today
            M = slice(P,j_mc,:)
            # *matrix* of continuous states tomorrow
            S = evaluate(g, m, s, x, M, p)
            # *matrix* of continous controls tomorrow
            X = evaluate(mdr, j_mc, S)
            res[i_mc,:,:] = slice(res,i_mc,:,:) + prob*evaluate(f,m,s,x,M,S,X,p)
            # res[i_mc,:,:] += evaluate(f,m,s,x,M,S,X,p)  # Dimension mismatch
        end

    end
    return res
end


function time_iteration(model::DTMSCCModel, mdrinit)

    f = model.functions.arbitrage
    g = model.functions.transition

    p = model.calibration[:parameters]

    # get approximation space
    ap = model.options[:approximation_space]
    a = ap[:a]
    b = ap[:b]
    dims = ap[:orders]
    approx = ApproximationSpace(a,b,dims)

    grid = mlinspace(a,b,dims)

    # get markov chain
    mcopts = model.options[:discrete_transition][:MarkovChain]
    nodes = cat(1, [e' for e  in mcopts[1]]...)
    transitions = cat(1, [e' for e  in mcopts[2]]...)
    markov_chain = MarkovChain(nodes, transitions)
    P = markov_chain.P
    Q = markov_chain.Q

    n_mc = size(P,1)
    n_m = size(P,2)

    n_x = length(model.calibration[:controls])
    N = size(grid,1)

    controls_0 = zeros(n_mc, N, n_x)
    for i_mc in 1:n_mc
        controls_0[i_mc,:,:] = evaluate(mdrinit, i_mc, grid)
    end
    s_x = size(controls_0)
    x0 = reshape(controls_0, n_mc*N, n_x)

    verbose = true

    mdr = MixedDecisionRule(n_mc,a,b,dims)

    # set iteration options
    maxit = 500
    tol = 1e-6

    err = 100
    it = 0

    while (err>tol)&&(it<maxit)
        it +=1
        set_values(mdr, reshape(x0, n_mc, N, n_x))
        fobj(t::Array{Float64,2}) = reshape( step_residuals(f, g, grid, reshape(t, s_x...), p, mdr, P, Q), size(x0)...)
        res = fobj(x0)
        (x,nit) = serial_solver(fobj, x0, 5)
        err = maximum(x-x0)
        x0 = x
        if verbose
            print((it, err, nit))
            print("\n")
        end
    end

    println("Finished in ", it, " iterations.")
    return mdr
end




####
# evaluate a given policy rule (a function)
####

function evaluate_policy(model::DTMSCCModel, policy, options=Dict())

    felicity = model.functions.felicity
    transition = model.functions.transition

    discount = model.calibration.flat[:beta]

    # get approximation space
    ap = model.options[:approximation_space]
    a = ap[:a]
    b = ap[:b]
    dims = ap[:orders]
    approx = ApproximationSpace(a,b,dims)

    grid = mlinspace(a,b,dims)
    d = size(grid,2)

    # get markov chain
    mcopts = model.options[:discrete_transition][:MarkovChain]
    nodes = cat(1, [e' for e  in mcopts[1]]...)
    transitions = cat(1, [e' for e  in mcopts[2]]...)
    markov_chain = MarkovChain(nodes, transitions)
    P = markov_chain.P
    Q = markov_chain.Q

    symbols = model.symbolic.symbols
    calibration = model.calibration

    N = size(grid,1)
    n_mc = size(markov_chain.P,1)
    n_s = size(grid,2)
    n_x = length(symbols[:controls])
    p = calibration
    controls = zeros(n_mc,N,n_x)
    value = zeros(n_mc,N)
    p = calibration[:parameters]

    for i_mc=1:n_mc
        for n=1:N
            s = slice(grid,n,:)
            m = slice(P,i_mc,:)
            x = evaluate(policy,i_mc,s)
            controls[i_mc, n,:] = x
            value[i_mc, n:n] *= evaluate(felicity, m, s, x, p)/(1-discount)
        end
    end

    it = 0
    maxit = get(options, :maxit, 500)
    tol =   get(options, :tol, 1e-6)
    err = 10

    dims = vcat(approx.dims)
    knots = [linspace(approx.a[i], approx.b[i], approx.dims[i]) for i=1:d]


    create_interpolant(vv) = [scale(interpolate(copy(reshape(slice(vv,i_mc,:),dims...)), BSpline(Linear()), OnGrid()), knots...)  for i_mc in 1:n_mc ]

    while it<maxit && err>tol

        it += 1
        value_0 = copy(value)


        fut_dr = create_interpolant(value_0)


        for i_mc=1:n_mc
            for n=1:N
                s = slice(grid, n, :)
                m = slice(P, i_mc, :)
                x = evaluate(policy,i_mc,s)
                val = evaluate(felicity,m, s, x, p)
                for j_mc in 1:n_mc
                    prob = Q[i_mc, j_mc]
                    # M = nodes[j_mc,:]
                    M = slice(P, j_mc, :)
                    S = evaluate(transition,m,s,x,M,p)
                    v_fut = fut_dr[j_mc][S...]
                    val += discount*prob*v_fut
                end
                value[i_mc, n:n] = val
            end
        end
        err = maximum( abs(value - value_0)[:] )
    end
    return value
end


####
# optimize policy rules by iterating on the value function
####

function solve_policy(model::DTMSCCModel, policy, values_0, options=Dict())

    felicity = model.functions.felicity
    transition = model.functions.transition

    discount = model.calibration.flat[:beta]

    # get approximation space
    ap = model.options[:approximation_space]
    a = ap[:a]
    b = ap[:b]
    dims = ap[:orders]
    approx = ApproximationSpace(a,b,dims)

    grid = mlinspace(a,b,dims)
    d = size(grid,2)

    # get markov chain
    mcopts = model.options[:discrete_transition][:MarkovChain]
    nodes = cat(1, [e' for e  in mcopts[1]]...)
    transitions = cat(1, [e' for e  in mcopts[2]]...)
    markov_chain = MarkovChain(nodes, transitions)
    P = markov_chain.P
    Q = markov_chain.Q

    symbols = model.symbolic.symbols
    calibration = model.calibration

    N = size(grid,1)
    n_mc = size(markov_chain.P,1)
    n_s = size(grid,2)
    n_x = length(symbols[:controls])
    p = calibration
    controls = zeros(n_mc,N,n_x)
    value = zeros(n_mc,N)
    p = calibration[:parameters]

    controls = zeros(n_mc,N,n_x)
    for i_mc=1:n_mc
        for n=1:N
            s = grid[n,:]
            m = P[i_mc,:]
            x = evaluate(policy,i_mc,s)
            controls[i_mc, n,:] = x
        end
    end

    lower_bound = model.functions.controls_lb
    upper_bound = model.functions.controls_ub

#     controls = copy(controls)
    value = copy(values_0)

    it = 0
    maxit = get(options, :maxit, 1000)
    tol = get(options, :tol, 0.00001)
    tol_2 = get(options, :tol_2, 0.00001)

    diff = 100
    diff_0 = 1.0
    diff_2_0 = 1.0

    dims = vcat(approx.dims)
    knots = [linspace(approx.a[i], approx.b[i], approx.dims[i]) for i=1:d]

    create_interpolant(vv) = [scale(interpolate(copy(reshape(slice(vv,i_mc,:),dims...)), BSpline(Cubic(Flat())), OnGrid()), knots...)  for i_mc in 1:n_mc ]


    while (it<maxit) & (diff>tol || diff_2>tol_2)

        tic()

        it += 1

        value_0 = copy(value)
        controls_0 = copy(controls)

        dims = vcat(approx.dims)

        fut_dr = create_interpolant(value_0)

        for i_mc=1:n_mc
            for n=1:N
                s = slice(grid,n,:)
                m = slice(P,i_mc,:)
                function fobj(x)
                    v = evaluate(felicity, m, s, x, p)
                    for j_mc in 1:n_mc
                        prob = slice(Q,i_mc, j_mc)
                        M = slice(P,j_mc,:)
                        S = evaluate(transition,m,s,x,M,p)

                        v_fut = fut_dr[j_mc][S...]

                        v += discount*prob*v_fut[1]
                    end
                    return -v[1]
                end

                x0 = copy(slice(controls, i_mc, n ,:))

                d4 = DifferentiableFunction(fobj)
                lb = evaluate(lower_bound,m,s,p)
                ub = evaluate(upper_bound,m,s,p)

                x0 = max( (min(x0, ub)), lb)

                # fminbox has apparently not been updated to the new Optim api
                res = fminbox(d4, x0, lb, ub)
                # res = optimize(d4, x0)

                value[i_mc, n:n] = -fobj(res.minimum)
                controls[i_mc, n,:] = res.minimum
            end
        end
        diff = maximum( abs(value - value_0)[:] )
        diff_2 = maximum( abs(controls - controls_0)[:])
        SA_r = diff/diff_0
        SA_r_2 = diff_2/diff_2_0
        diff_0 = diff
        diff_2_0 = diff_2
        elapsed = toq()
        println("(", it,", ",diff,", ", SA_r, ", ", diff_2, ", ", SA_r_2, " : ", elapsed, ")")
    end
    return controls,value
end
