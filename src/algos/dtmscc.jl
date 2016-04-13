using QuantEcon: gridmake

using Optim

using Interpolations
using splines # temporary

#####
# helper functions/types (to be moved/removed)
#####

mlinspace(a,b,orders) = gridmake([linspace(a[i],b[i],orders[i]) for i=1:length(orders)]...)


type MarkovChain
    P::Array{Float64, 2}
    Q::Array{Float64, 2}
end

type ApproximationSpace
  a:: Array{Float64, 1}
  b:: Array{Float64, 1}
  orders:: Array{Int64, 1}
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

function constant_guess(model::DTMSCCModel)
    x = model.calibration[:controls]
    function fun(i,s)
        if ndims(s) == 1
            return x
        else
            N = size(s,1)
            return repmat(x',N,1)
        end
    end
    return fun
end

####
# evaluate a given policy rule (a function)
####

function eval_policy(model::DTMSCCModel, policy, options=Dict())

    felicity = model.functions.felicity
    transition = model.functions.transition

    discount = model.calibration.flat[:beta]

    # get approximation space
    ap = model.options[:approximation_space]
    a = ap[:a]
    b = ap[:b]
    orders = ap[:orders]
    approx = ApproximationSpace(a,b,orders)

    grid = mlinspace(a,b,orders)
    d = size(grid,2)

    # get markov chain
    mcopts = model.options[:discrete_transition][:MarkovChain]
    P = cat(1, [e' for e  in mcopts[1]]...)
    Q = cat(1, [e' for e  in mcopts[2]]...)
    markov_chain = MarkovChain(P,Q)
    nodes = markov_chain.P

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
            m = slice(nodes,i_mc,:)
            x = policy(i_mc,s)
            controls[i_mc, n,:] = x
            value[i_mc, n:n] *= evaluate(felicity, m, s, x, p)/(1-discount)
        end
    end

    print(value)
    it = 0
    maxit = get(options, :maxit, 500)
    tol =   get(options, :tol, 1e-6)
    err = 10

    dims = vcat(approx.orders)
    knots = [collect(linspace(approx.a[i], approx.b[i], approx.orders[i])) for i=1:d]
    knots = tuple(knots...)

    create_interpolant(vv) = [interpolate(knots, copy(reshape(slice(vv,i_mc,:),dims...)), Gridded(Linear()))  for i_mc in 1:n_mc ]

    while it<maxit && err>tol

        it += 1
        value_0 = copy(value)


        fut_dr = create_interpolant(value_0)


        for i_mc=1:n_mc
            for n=1:N
                s = slice(grid, n, :)
                m = slice(nodes, i_mc, :)
                x = policy(i_mc,s)
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
    orders = ap[:orders]
    approx = ApproximationSpace(a,b,orders)

    grid = mlinspace(a,b,orders)
    d = size(grid,2)

    # get markov chain
    mcopts = model.options[:discrete_transition][:MarkovChain]
    P = cat(1, [e' for e  in mcopts[1]]...)
    Q = cat(1, [e' for e  in mcopts[2]]...)
    markov_chain = MarkovChain(P,Q)
    nodes = markov_chain.P

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
            m = nodes[i_mc,:]
            x = policy(i_mc,s)
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

    dims = vcat(approx.orders)
    knots = [collect(linspace(approx.a[i], approx.b[i], approx.orders[i])) for i=1:d]
    knots = tuple(knots...)

    # interpolations.jl does not have cubic gridded interpolations yet
    # using splnes.jl isntead in the meantime
     create_interpolant(vv) = [interpolant_cspline(approx.a, approx.b, approx.orders, copy(reshape(slice(vv,i_mc,:),dims...) )) for i_mc in 1:n_mc]
    # create_interpolant(vv) = [ interpolate(knots, copy(reshape(slice(vv,i_mc,:),dims...)), Gridded(Cubic()))  for i_mc in 1:n_mc ]

    while (it<maxit) & (diff>tol || diff_2>tol_2)

        tic()

        it += 1

        value_0 = copy(value)
        controls_0 = copy(controls)

        dims = vcat(approx.orders)

        fut_dr = create_interpolant(value_0)

        for i_mc=1:n_mc
            for n=1:N
                s = slice(grid,n,:)
                m = slice(nodes,i_mc,:)
                function fobj(x)
                    v = evaluate(felicity, m, s, x, p)
                    for j_mc in 1:n_mc
                        prob = slice(Q,i_mc, j_mc)
                        M = slice(nodes,j_mc,:)
                        S = evaluate(transition,m,s,x,M,p)

                        # v_fut = fut_dr[j_mc][S...]
                        v_fut = fut_dr[j_mc](S...)

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
