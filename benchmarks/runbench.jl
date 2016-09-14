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

const git_dir = joinpath(dirname(dirname(@__FILE__)), ".git")
const sha = readchomp(`git --git-dir $(git_dir) rev-parse HEAD`)

if !isfile("results.jld")
    close(jldopen("results.jld", "w"))
end

jldopen("results.jld", "r+") do f
    nm = "results_$(sha)"
    JLD.exists(f, nm) && JLD.delete!(f, nm)
    JLD.write(f, nm, results)
end

end  # module
