using TimerOutputs: TimerOutput


struct IterationLog
    header::Array{String, 1}
    keywords::Array{Symbol, 1}
    types::Vector
    entries::Array{Any, 1}
    fmt::String
end

function IterationLog(;kw...)
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
    IterationLog(header, keywords, types, [], fmt)
end

function append!(log::IterationLog; verbose=true, entry...)
    d = Dict(entry)
    push!(log.entries, d)
    verbose && show_entry(log, d)
end

function initialize(log::IterationLog; verbose=true, message=nothing)
    verbose && show_start(log; message=message)
end

function finalize(log::IterationLog; verbose=true)
    verbose && show_end(log)
end

function show(log::IterationLog)
    show_start(log)
    for entry in log.entries
        show_entry(log,entry)
    end
    show_end(log)
end

function header(log::IterationLog)
    l = []
    for t in log.types
        if t<:Int
            push!(l, "{:>8}")
        elseif t<:Real
            push!(l, "{:>12}")
        end
    end
    return format(join(l, " | "), log.header...)
end

function show_start(log::IterationLog; message=nothing)

    s = header(log)

    if message !== nothing
        println(repeat("-", length(s)))
        println(message)
    end
    println(repeat("-", length(s)))

    println(s)
    println(repeat("-", length(s)))
end

function show_entry(log::IterationLog, entry)

    vals = [entry[k] for k in log.keywords]
    printfmt(log.fmt, vals...)
    print("\n")

end

function show_end(log::IterationLog)
    println(repeat("-", length(header(log))))
end

###


struct IterationTrace
    data::Array{Any,1}
end

length(t::IterationTrace) = length(t.data)


struct TimeIterationResult
    dr
    iterations::Int
    tol_ε::Float64
    ε::Float64
    tol_η::Float64
    η::Float64
    log::Union{Nothing, IterationLog}
    trace::Union{Nothing,IterationTrace}
    timer::Union{Nothing,TimerOutput}
end

converged(r::TimeIterationResult) = (r.ε<r.tol_ε) || (r.η<r.tol_η)
function Base.show(io::IO, r::TimeIterationResult)
    @printf io "Results of Time Iteration Algorithm\n"
    @printf io " * Decision Rule type: %s\n" string(typeof(r.dr))
    @printf io " * Number of iterations: %s\n" string(r.iterations)
    @printf io " * Convergence: %s\n" converged(r) ? "true" : RED_FG("false")
    @printf io "   * |F(x,x)| < %.1e: %s (%.1e)\n" r.tol_ε (r.ε<r.tol_ε) r.ε
    @printf io "   * |x - x'| < %.1e: %s (%.1e)\n" r.tol_η (r.η<r.tol_η) r.η
end


struct NewtonResult
    dr
    iterations::Int
    tol_ε::Float64
    ε::Float64
    log::Union{Nothing, IterationLog}
    timer::Union{Nothing,TimerOutput}
end

converged(r::NewtonResult) = (r.ε<r.tol_ε)
function Base.show(io::IO, r::NewtonResult)
    @printf io "Results of Newton Algorithm\n"
    @printf io " * Decision Rule type: %s\n" string(typeof(r.dr))
    @printf io " * Number of iterations: %s\n" string(r.iterations)
    @printf io " * Convergence: %s\n" converged(r) ? "true" : RED_FG("false")
    @printf io "   * |F(x,x)| < %.1e: %s (%.1e)\n" r.tol_ε converged(r) r.ε
end


######
###### Value iteration
######


struct ValueIterationLog
    header::Array{String, 1}
    keywords::Array{Symbol, 1}
    entries::Array{Any, 1}
end

function ValueIterationLog()
    header = [ "It", "ηₙ=|xₙ-xₙ₋₁|", "νₙ=|vₙ-vₙ₋₁|", "Time", "Eval steps"]
    keywords = [:it, :sa_x, :sa_v, :time, :nit]
    ValueIterationLog(header, keywords, [])
end


struct ValueIterationResult
    dr
    drv
    iterations::Int
    x_tol::Float64
    x_err::Float64
    v_tol::Float64
    v_err::Float64
    log::ValueIterationLog
    trace::Union{Nothing,IterationTrace}
end

converged(r::ValueIterationResult) = r.x_err < r.x_tol
function Base.show(io::IO, r::ValueIterationResult)
    @printf io "Results of Value Iteration Algorithm\n"
    @printf io " * Decision Rule type: %s\n" string(typeof(r.dr))
    @printf io " * Number of iterations: %s\n" string(r.iterations)
    @printf io " * Convergence: %s\n" converged(r)
    @printf io "   * |x - x'| < %.1e: %s\n" r.x_tol converged(r)
    # @printf io "   * |v - v'| < %.1e: %s\n" r.v_tol r.v_converged

end

function append!(log::ValueIterationLog; verbose=true, entry...)
    d = Dict(entry)
    push!(log.entries, d)
    verbose && show_entry(log, d)
end

function initialize(log::ValueIterationLog; verbose=true)
    verbose && show_start(log)
end

function finalize(log::ValueIterationLog; verbose=true)
    verbose && show_end(log)
end

function show(log::ValueIterationLog)
    show_start(log)
    for entry in log.entries
        show_entry(log,entry)
    end
    show_end(log)
end

function show_start(log::ValueIterationLog)
    println(repeat("-", 66))
    @printf "%-6s%-16s%-16s%-16s%-5s\n" "It" "|xₙ-xₙ₋₁|" "|vₙ-vₙ₋₁|" "Time" "Eval steps"
    println(repeat("-", 66))
end


function show_entry(log::ValueIterationLog, entry)
    it = entry[:it]
    sa_x = entry[:sa_x]
    sa_v = entry[:sa_v]
    time = entry[:time]
    nit = entry[:nit]
    @printf "%-6i%-16.2e%-16.2e%-16.2e%-5i\n" it sa_x sa_v time nit
end

function show_end(log::ValueIterationLog)
    println(repeat("-", 66))
end
