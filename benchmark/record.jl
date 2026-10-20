# Results of the benchmarks, saved with the platform, the software and the commit, so that
# runs can be compared across machines and over time (see compare.jl).
#
# One TOML file per run, in benchmark/history/ (tracked by git), named
# `<date>_<host>_<script>_<backend>.toml`:
#
#     [run]        script, date, options
#     [platform]   host, os, cpu, threads, memory, gpu
#     [software]   julia, Dolo commit, versions of the packages, GPU driver
#     [results]    rows = [ one line per measurement: benchmark, model, precision, N, grid,
#                  engine, operation, time (seconds), status, and optional fields (iterations,
#                  residual, maxdiff, error, sections) ]
#
# Read them with `TOML.parsefile` (the measurements are in `d["results"]["rows"]`), or with
# `load_history` below.

using TOML, LibGit2, Dates, Printf

const HISTORY = joinpath(@__DIR__, "history")

"Commit of the Dolo.jl checkout in `dir` (read with LibGit2: no `git` executable needed)."
function git_info(dir)
    try
        repo = LibGit2.GitRepo(dir)
        oid = LibGit2.head_oid(repo)
        commit = LibGit2.GitCommit(repo, oid)
        Dict(
            "commit" => string(oid),
            "branch" => LibGit2.headname(repo),
            "dirty" => LibGit2.isdirty(repo),
            "subject" => first(split(LibGit2.message(commit), '\n')),
            "date" => string(unix2datetime(LibGit2.author(commit).time)),
        )
    catch e
        Dict("commit" => "unknown", "error" => sprint(showerror, e))
    end
end

function package_versions(names)
    d = Dict{String,String}()
    for (id, m) in Base.loaded_modules
        id.name in names || continue
        v = pkgversion(m)
        v === nothing || (d[id.name] = string(v))
    end
    return d
end

kernel_release() = Sys.isunix() ? strip(read(`uname -r`, String)) : ""

"Metadata of a run: script, options, platform, software and commit."
function run_metadata(script; gpu=nothing, options=Dict{String,Any}())
    platform = Dict{String,Any}(
        "host" => gethostname(),
        "os" => "$(Sys.KERNEL) $(kernel_release()) $(Sys.MACHINE)",
        "cpu" => Sys.cpu_info()[1].model,
        "cpu_threads" => Sys.CPU_THREADS,
        "julia_threads" => Threads.nthreads(),
        "memory_GiB" => round(Sys.total_memory() / 2^30, digits=1),
    )
    software = Dict{String,Any}(
        "julia" => string(VERSION),
        "dolo" => git_info(pkgdir(Dolo)),
        "packages" => package_versions(["Dolo", "KernelAbstractions", "GPUArrays", "GPUCompiler",
            "StaticArrays", "CUDA", "cuBLAS", "oneAPI", "LLVM"]),
    )
    if gpu !== nothing && gpu.functional
        platform["gpu"] = gpu.platform
        software["gpu"] = gpu.software
    end
    run = Dict{String,Any}("script" => script, "date" => string(now()), "options" => options)
    backend = gpu === nothing ? "cpu" : gpu.name
    path = joinpath(HISTORY, "$(Dates.format(now(), "yyyy-mm-dd_HHMM"))_$(platform["host"])_$(script)_$(backend).toml")
    return (; path, meta=Dict{String,Any}("run" => run, "platform" => platform, "software" => software))
end


## Writing (the results are written by hand: one inline table per line)

toml_string(s) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\r" => "",
    "\t" => "\\t", r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]" => "") * "\""
toml_value(x::AbstractString) = toml_string(x)
toml_value(x::Symbol) = toml_string(string(x))
toml_value(x::Bool) = string(x)
toml_value(x::Integer) = string(x)
toml_value(x::AbstractFloat) = isnan(x) ? "nan" : isinf(x) ? (x > 0 ? "inf" : "-inf") : repr(Float64(x))
toml_value(x::AbstractDict) = "{" * join(("$(toml_string(k)) = $(toml_value(v))" for (k, v) in x), ", ") * "}"
toml_value(x::NamedTuple) = "{" * join(("$k = $(toml_value(v))" for (k, v) in pairs(x) if v !== nothing), ", ") * "}"
toml_value(x::DataType) = toml_string(string(x))

"Write the metadata and the results (a vector of named tuples) of a run. Called after each step, so
that nothing is lost if a run is interrupted."
function save(run, rows)
    mkpath(dirname(run.path))
    open(run.path, "w") do io
        TOML.print(io, run.meta; sorted=true)
        println(io, "\n[results]\nrows = [")
        for r in rows
            println(io, "    ", toml_value(r), ",")
        end
        println(io, "]")
    end
    return run.path
end

"Short description of an error, for the tables."
function short_error(e)
    msg = sprint(showerror, e)
    m = match(r"Total arguments size: (\d+), limit: (\d+)", msg)
    m !== nothing && return "kernel arguments $(m[1]) B > $(m[2]) B"
    occursin("Out of GPU memory", msg) && return "out of GPU memory"
    # oneAPI: the cause (e.g. kernel arguments too large) is in the build log, printed before
    occursin("ZE_RESULT_ERROR_MODULE_BUILD_FAILURE", msg) && return "kernel compilation failed"
    return first(replace(msg, r"\s+" => " "), 200)
end

"Times of the sections of a TimerOutput, in seconds, with the nested ones as `outer / inner`."
function timer_sections(tmr, prefix="")
    d = Dict{String,Float64}()
    tmr === nothing && return d
    for (name, t) in tmr.inner_timers
        key = prefix * name
        d[key] = t.accumulated_data.time / 1e9
        merge!(d, timer_sections(t, key * " / "))
    end
    return d
end


## Reading and tables

"All the results in the TOML files `paths`, each row with the metadata of its run."
function load_history(paths=filter(endswith(".toml"), readdir(HISTORY; join=true)))
    rows = Dict{String,Any}[]
    for path in paths
        d = TOML.parsefile(path)
        host = d["platform"]["host"]
        commit = get(d["software"]["dolo"], "commit", "unknown")
        dirty = get(d["software"]["dolo"], "dirty", false)
        for r in get(get(d, "results", Dict()), "rows", [])
            push!(rows, merge(r, Dict("host" => host, "date" => d["run"]["date"],
                "commit" => first(commit, 7) * (dirty ? "+" : ""), "file" => basename(path), "path" => path)))
        end
    end
    return rows
end

format_time(t) = !isfinite(t) ? "" : t < 1e-3 ? @sprintf("%.0f μs", t*1e6) : t < 1 ? @sprintf("%.1f ms", t*1e3) : @sprintf("%.2f s", t)

const ENGINE_ORDER = ["serial", "threads", "cuda", "oneapi"]
const PRECISION_ORDER = ["Float64", "Float32"]
const OPERATION_ORDER = ["fit!", "F!", "dF_1!", "dF_2!", "L.v", "J\\L", "time iteration", "accelerated TI", "newton"]

function cell(r)
    get(r, "status", "ok") == "error" && return "✗ " * get(r, "error", "error")
    s = format_time(r["time"])
    haskey(r, "iterations") && (s *= " ($(r["iterations"]) it)")
    get(r, "status", "ok") == "no convergence" && (s *= " ✗ no conv.")
    return s
end

"""
Markdown table of `rows` (dictionaries with string keys): one line per benchmark, precision, grid
and operation, one column per `column(row)` (by default host and engine). When several rows fall
in the same cell, the latest run wins.
"""
function print_table(io::IO, rows; column=r -> "$(r["host"]) $(r["engine"])")
    isempty(rows) && return
    rows = sort(rows; by=r -> r["date"])
    cols = sort(unique(column.(rows)); by=c -> (split(c)[1], something(findfirst(==(split(c)[end]), ENGINE_ORDER), 99), c))
    rowkey(r) = (r["benchmark"], r["precision"], r["N"], r["operation"])
    keys_ = sort(unique(rowkey.(rows)); by=k -> (k[1], something(findfirst(==(k[2]), PRECISION_ORDER), 99), k[2], k[3],
        something(findfirst(==(k[4]), OPERATION_ORDER), 99)))
    cells = Dict((rowkey(r), column(r)) => cell(r) for r in rows)
    grids = Dict(r["N"] => r["grid"] for r in rows)
    println(io, "| benchmark | precision | grid | operation | ", join(cols, " | "), " |")
    println(io, "|---|---|---|---|", repeat("---|", length(cols)))
    for k in keys_
        println(io, "| ", k[1], " | ", k[2], " | ", grids[k[3]], " | ", k[4], " | ",
            join((get(cells, (k, c), "") for c in cols), " | "), " |")
    end
end
print_table(rows; kw...) = print_table(stdout, rows; kw...)

"Rows of the current run (named tuples) in the form of `load_history`."
function as_history(run, rows)
    [merge(Dict(string(k) => v for (k, v) in pairs(r) if v !== nothing),
        Dict("host" => run.meta["platform"]["host"], "date" => run.meta["run"]["date"]))
     for r in rows]
end
