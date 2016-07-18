module DoloBenchmarks
using Dolo
using Dolang
using BenchmarkTools
using JLD

suite = BenchmarkGroup()

include("ea_quest.jl")


if "tune" in ARGS || !isfile("params.jld")
    println("Tuning benchmarks")
    tune!(suite)
    JLD.save("params.jld", "suite", params(suite))
else
    println("loading saved tuning params")
    loadparams!(suite, JLD.load("params.jld", "suite"), :evals, :samples)
end

results = run(suite, verbose=true)

JLD.save("results.jld", "results", results)

end  # module
