
struct Iterationlog
    header::Array{String, 1}
    keywords::Array{Symbol, 1}
    types::Vector
    entries::Array{Any, 1}
    fmt::String
end

# function Iterationlog()
#     header = [ "It", "ϵₙ", "ηₙ=|xₙ-xₙ₋₁|", "λₙ=ηₙ/ηₙ₋₁", "Time", "Newton steps"]
#     keywords = [:it, :err, :gain, :time, :nit]
#     Iterationlog(header, keywords, [])
# end

function Iterationlog(;kw...)
    keywords = Symbol[]
    header = String[]
    elts = []
    types= []
    for (k,(v,t)) in kw
        push!(types, t)
        push!(keywords, k)
        push!(header, v)
        if t<:Int
            push!(elts, "{:>8}")
        elseif t<:Real
            push!(elts, "{:>12.4e}")
        else
            error("Don't know how to format type $t")
        end
    end
    fmt = join(elts, " | ")
    Iterationlog(header, keywords, types, [], fmt)
end

function append!(log::Iterationlog; verbose=true, entry...)
    d = Dict(entry)
    push!(log.entries, d)
    verbose && show_entry(log, d)
end

function initialize(log::Iterationlog; verbose=true, message=nothing)
    verbose && show_start(log; message=message)
end

function finalize(log::Iterationlog; verbose=true)
    verbose && show_end(log)
end

function show(log::Iterationlog)
    show_start(log)
    for entry in log.entries
        show_entry(log,entry)
    end
    show_end(log)
end

function show_start(log::Iterationlog; message=nothing)
    # println(repeat("-", 66))
    # @printf "%-6s%-16s%-16s%-16s%-16s%-5s\n" "It" "ϵₙ" "ηₙ=|xₙ-xₙ₋₁|" "λₙ=ηₙ/ηₙ₋₁" "Time" "Newton steps"
    # println(repeat("-", 66))

    l = []
    for t in log.types
        if t<:Int
            push!(l, "{:>8}")
        elseif t<:Real
            push!(l, "{:>12}")
        end
    end

    a = join(l, " | ")
    s = format(a, log.header...)

    if message !== nothing
        println(repeat("-", length(s)))
        println(message)
    end
    println(repeat("-", length(s)))

    println(s)
    println(repeat("-", length(s)))
end


function show_entry(log::Iterationlog, entry)

    vals = [entry[k] for k in log.keywords]
    printfmt(log.fmt, vals...)
    print("\n")

end

function show_end(log::Iterationlog)
    println(repeat("-", 66))
end

###

mutable struct IterationTrace
    trace::Array{Any,1}
end

mutable struct TimeIterationResult
    dr::AbstractDecisionRule
    iterations::Int
    complementarities::Bool
    dprocess::AbstractDiscretizedProcess
    ϵ::Float64 
    η::Float64
    τ_ϵ::Float64 
    τ_η::Float64
    λ_bar::Float64 
    p::Int 
    log::Iterationlog
    trace::Union{Nothing,IterationTrace}
end


solved(r::TimeIterationResult) = r.ϵ<r.τ_ϵ       
eta_converged(r::TimeIterationResult) = r.η<r.τ_η        


function Base.show(io::IO, r::TimeIterationResult)
    @printf io "Results of Time Iteration Algorithm\n"
    @printf io " * Complementarities: %s\n" string(r.complementarities)
    @printf io " * Discretized Process type: %s\n" string(typeof(r.dprocess))
    @printf io " * Decision Rule type: %s\n" string(typeof(r.dr))
    @printf io " * Number of iterations: %s\n" string(r.iterations)
    @printf io " * ϵₙ=|F(xₙ,xₙ)| < %.1e: %s\n" r.τ_ϵ solved(r)
    @printf io " * ηₙ=|x - x'| < %.1e: %s\n" r.τ_η eta_converged(r)
end


function newton2(fun, dfun, x0::MSM; maxit=50, dampen=1.0, verbose=false, bcksteps = 5)

    r0 = fun(x0)

    err = norm(r0)

    local x1

    x1data = deepcopy(x0.data)
    x1 = MSM(x1data, x0.sizes)


    for n=1:maxit

        r0 = fun(x0)
        err_0 = norm(r0)

        j = dfun(x0)
        δ = (j\r0)
        for i=0:(bcksteps-1)
            u = 0.5^i
            x1.data[:] .= x0.data - δ.data*u
            r1 = fun(x1)
            err_1 = norm(r1)
            if verbose
                println( "-    $i: ", err_1)
            end
            if err_1<err_0
                break
            end
        end

        if verbose
            println(n, " | ", norm(r0), " | ",  norm(δ))
        end

        x0 = x1

    end

    return x0

end


"""
Computes a global solution for a model via backward time iteration. The time iteration is applied to the residuals of the arbitrage equations.

If the initial guess for the decision rule is not explicitly provided, the initial guess is provided by `ConstantDecisionRule`.
If the stochastic process for the model is not explicitly provided, the process is taken from the default provided by the model object, `model.exogenous`

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `process`: The stochastic process associated with the exogenous variables in the model.
* `init_dr`: Initial guess for the decision rule.
# Returns
* `dr`: Solved decision rule.
"""


function time_iteration(model;
    dr0=Dolo.ConstantDecisionRule(model),
    discretization=Dict(),
    interpolation=:cubic,
    verbose=true,
    details=true,
    ignore_constraints=false,
    trace = false,
    tol_η = 1e-8,
    tol_ε = 1e-8,
    λbar = 8.0466e-01,   
    maxit=500
)
    
    F = Euler(model; discretization=discretization, interpolation=interpolation, dr0=dr0,  ignore_constraints=ignore_constraints)

    complementarities = false

    if interpolation != :cubic
        error("Interpolation option ($interpolation) is currently not recognized.")
    end

    z0 = deepcopy(F.x0)
    z1 = deepcopy(z0)

    ti_trace = trace ? IterationTrace([deepcopy(F.dr.dr)]) : nothing

    vector_time = rand(maxit) 

    local err_ε, err_η, z0, z1, p, it 

    log = Iterationlog(
        it = ("n", Int),
        err =  ("εₙ=|F(xₙ,xₙ)|", Float64),
        sa =  ("ηₙ=|xₙ-xₙ₋₁|", Float64),
        lam = ("λₙ=ηₙ/ηₙ₋₁", Float64),
        elapsed = ("Time (s)", Float64),
        nb_it_before_convergence_of_x = ("N-n", Int), 
        remaining_time = ("Remain (s)", Float64), 
    )

    initialize(log, verbose=verbose; message="Time Iteration")

    err_η_0 = NaN

    it = 0

    while it<=maxit

        it += 1

        t1 = time_ns()

        r0 = F(z0, z0; set_future=true)

        err_ε=norm(r0)

        fun = u->F(u, z0; set_future=false)
        dfun = u->df_A(F,u, z0; set_future=false)

        sol = newton2(
            fun, dfun,
            z0, verbose=false
        )

        z1 = sol
        δ = z1 - z0

        trace && push!(ti_trace.trace, deepcopy(F.dr.dr)) 

        err_η = norm(δ)
        gain = err_η / err_η_0
        err_η_0 = err_η

        # z0.data[:] .= z1.data
        z0 = z1
        
        p = NaN

        if gain>=1
            p = NaN
        elseif gain<1
            p = Base.log(tol_η/err_η) / Base.log(λbar)
        end

        elapsed = time_ns() - t1

        elapsed /= 1000000000

        vector_time[it] = elapsed 

        local avg_time
        
        avg_time = mean(vector_time[max(1,end-4):end])

        time_left = NaN
        if gain>=1
            time_left = NaN
            p = NaN
        elseif gain<1
            time_left = avg_time*p
            p = round(Int, p) 
        end

     append!(log; 
        verbose=verbose, 
        it=it, 
        err=err_ε, 
        sa=err_η, 
        lam=gain, 
        elapsed=elapsed,
        nb_it_before_convergence_of_x=p, 
        remaining_time=time_left,
     )

     if err_ε<tol_ε     
         break
     end

     if err_η<tol_η       
         break 
     end

end

 finalize(log, verbose=verbose)

 res = TimeIterationResult(
     F.dr.dr, 
     it, 
     complementarities, 
     F.dprocess, 
     err_ε, 
     err_η, 
     tol_ε, 
     tol_η, 
     λbar, 
     p, 
     log, 
     ti_trace
     )

 return res

end
 

function improved_time_iteration(model;
    dr0=Dolo.ConstantDecisionRule(model),
    discretization=Dict(),
    interpolation=:cubic,
    verbose=true,
    details=true,
    ignore_constraints=false,
    trace = false,
    tol_η = 1e-8,
    tol_ε = 1e-8,
    tol_ν = 1e-10,
    maxit=500,
    smaxit=500
)


    F = Euler(model; discretization=discretization, interpolation=interpolation, dr0=dr0, ignore_constraints=ignore_constraints)

    complementarities = false


    if interpolation != :cubic
        error("Interpolation option ($interpolation) is currently not recognized.")
    end

    ti_trace = trace ? IterationTrace([x0]) : nothing


    z0 = deepcopy(F.x0)

    local err_η, err_ε

    log = IterationLog(;
        it=("It",Int),
        err= ("ϵₙ=|F(xₙ,xₙ)|",Float64),
        sa= ("ηₙ=|xₙ-xₙ₋₁|", Float64),
        time=("Time", Float64)
    )
    initialize(log; verbose=verbose, message="Improved Time Iteration")

 
    it = 0
    while it<=maxit

        t1 = time_ns()

        it += 1

        r = F(z0, z0; set_future=true)

        err_ε = norm(r)

        J = df_A(F, z0, z0 ; set_future=false)
        L = df_B(F, z0, z0 ; set_future=false)

        mult!(L, -1.0)
        prediv!(L, J) # J\L

        π = -J\r

        count = 0

        u = π
        δ = π
        for i=1:smaxit
            count +=1
            u = L*u
            δ += u
            if norm(u)<tol_ν
                break
            end
        end
        

        z0 = z0 + δ

        err_η = norm(δ)

        elapsed = time_ns() - t1

        append!(log; verbose=verbose,
            it = it, 
            err = err_ε,
            sa = err_η,
            time = elapsed
        )

        if err_η<tol_η
            break 
        end

    end

    finalize(log, verbose=verbose)

    res = ImprovedTimeIterationResult(
        F.dr.dr,
        it,
        err_ε,
        err_η,
        err_η<tol_η,
        complementarities,
        F.dprocess,
        tol_ν,
        NaN,
        0,
        0,
        0,
        ti_trace
    )
        
    res
end



function update_guess(F, x0, p0, dp, tol=1e-8, maxit=1000)
    df = dFdp(F, x0, x0, p0, dp)
    J = Dolo.df_A(F, x0, x0)
    L = Dolo.df_B(F, x0, x0)
    Dolo.prediv!(L, J) 
    Dolo.mult!(L, -1.0)
    dπ = -J\df
    dx = dπ
    i = 0
    while i<=maxit
        i += 1
        err = Dolo.norm(dπ)
        if err<tol
            break
        end
        dπ = L*dπ
        dx += dπ
    end
    return dx, i
end
