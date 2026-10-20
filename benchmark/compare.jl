# Compare the results saved in benchmark/history/: one table, one column per machine and engine.
#
#     julia --project=benchmark benchmark/compare.jl [files...] [--benchmark=kernels|solve]
#                    [--precision=Float32] [--grid=20³×9] [--hosts=takeshi,parker] [--by=commit]
#
# By default all the files of benchmark/history/ are read and, for each cell, the latest run wins.
# With `--by=commit`, each commit gets its own columns (for comparisons over time).

include("record.jl")

option(name) = let a = findfirst(startswith("--$name="), ARGS)
    a === nothing ? nothing : split(split(ARGS[a], "=")[2], ",")
end

files = filter(a -> !startswith(a, "--"), ARGS)
rows = load_history(isempty(files) ? filter(endswith(".toml"), readdir(HISTORY; join=true)) : files)

for (key, field) in [("benchmark", "benchmark"), ("precision", "precision"), ("grid", "grid"), ("hosts", "host")]
    values = option(key)
    values === nothing || filter!(r -> r[field] in values, rows)
end

column = option("by") == ["commit"] ? (r -> "$(r["host"])@$(r["commit"]) $(r["engine"])") : (r -> "$(r["host"]) $(r["engine"])")

println("Runs:\n")
for path in sort(unique(r["path"] for r in rows); by=basename)
    file = basename(path)
    d = TOML.parsefile(path)
    p, s = d["platform"], d["software"]
    gpu = haskey(p, "gpu") ? ", $(p["gpu"]["name"]) ($(s["gpu"]["backend"]) $(get(s["gpu"], "driver", "")))" : ""
    dolo = s["dolo"]
    println("- ", file, ": ", p["cpu"], ", ", p["julia_threads"], " threads", gpu, "; Julia ", s["julia"],
        "; Dolo ", first(get(dolo, "commit", "unknown"), 7), get(dolo, "dirty", false) ? " (modified)" : "")
end
println()
print_table(rows; column)
