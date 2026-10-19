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
# complete solution: 5 time iteration steps, then Newton (GMRES), on threads and GPU
julia -t auto --project=benchmark benchmark/integration_B_solve.jl --sizes=20
```

Both take `--backend=cuda` (default) or `--backend=oneapi` (Intel GPUs), which must then be
run in the `benchmark/oneapi` environment:

```bash
julia --project=benchmark/oneapi -e 'using Pkg; Pkg.instantiate()'
julia -t auto --project=benchmark/oneapi benchmark/integration_B_solve.jl --backend=oneapi
```

On GPUs without Float64 (many Intel GPUs), the GPU runs in Float32 only.

