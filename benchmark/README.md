# Benchmarks

Performance of Dolo.jl across models, grid sizes and engines (serial, multithreaded, GPU).

```bash
julia -t auto --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia -t auto --project=benchmark benchmark/run.jl
```

Options:

- `--quick`: smaller grids, for a fast check
- `--models=rbc_mc,rbc_ar1`: models, from `examples/ymodels`
- `--engines=serial,threads,gpu`: engines (`gpu` only if CUDA is functional)
- `--only=kernels` or `--only=solvers`: one category of benchmarks

Two categories are measured:

- **kernels**: one call of each vectorized operation on the whole grid:
  `fit!` (spline fitting of the decision rule), `F!` (residuals), `dF_1!` and `dF_2!`
  (derivatives), `L.v` (application of the linear operator, the inner loop of
  Neumann and GMRES) and `J\L`. Grids of 10³ to 10⁶ points.
- **solvers**: complete solution from the initial guess, with plain time iteration,
  accelerated time iteration and Newton (with backtracking). Grids of 10³ and 10⁴ points.

All use cubic splines and Float64. Results are printed as markdown tables, with speed-ups
relative to the serial engine, and saved in `benchmark/results/` (markdown and CSV, not
tracked by git).

## integration_B

Two scripts focus on the integration_B model (8 controls, 9 exogenous nodes, 3 endogenous
states), which is large enough for the GPU to pay off:

```bash
# one call of each operation, on CPU (serial, threads) and GPU, in Float64 and Float32
julia -t auto --project=benchmark benchmark/integration_B.jl --sizes=10,20,30,40
# complete solution from the initial guess, on threads and GPU, with time iteration,
# accelerated time iteration and Newton (with backtracking)
julia -t auto --project=benchmark benchmark/integration_B_solve.jl --sizes=20
```

Both take `--backend=cuda` (default) or `--backend=oneapi` (Intel GPUs), which must then be
run in the `benchmark/oneapi` environment:

```bash
julia --project=benchmark/oneapi -e 'using Pkg; Pkg.instantiate()'
julia -t auto --project=benchmark/oneapi benchmark/integration_B_solve.jl --backend=oneapi --engines=gpu
```

On GPUs without Float64 (many Intel GPUs), the GPU runs in Float32 only. A kernel or an
algorithm that fails (out of memory, kernel arguments too large...) is recorded as an error,
and the others still run.

## History of the results

The integration_B scripts save their results in `benchmark/history/` (tracked by git), one
TOML file per run, `<date>_<host>_<script>_<backend>.toml`, with:

- `[platform]`: host, OS, CPU, threads, memory, GPU (name, memory, Float64 support)
- `[software]`: Julia, Dolo commit (hash, branch, modified or not), versions of the GPU
  packages, GPU driver and runtime
- `[results]`: one line per measurement (`benchmark`, `precision`, `grid`, `engine`,
  `operation`, `time` in seconds, `status`; `iterations`, `residual` and the times of the
  sections of the solvers for `solve`; `maxdiff` from the CPU for the GPU kernels; `error`)

To compare runs (machines, commits), in one table:

```bash
julia --project=benchmark benchmark/compare.jl                       # latest result of each machine and engine
julia --project=benchmark benchmark/compare.jl --benchmark=solve --precision=Float32
julia --project=benchmark benchmark/compare.jl --by=commit --hosts=parker   # over time
```

Commit before running the benchmarks, so that the results refer to a commit
(`modified` otherwise).
