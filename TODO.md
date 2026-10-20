# TODO

## oneAPI (Intel GPU, on the laptop)

Prepared on the `integration-gpu` branch, but not run on an Intel GPU yet: the desktop's
UHD 630 (Gen9) is not supported by the Intel driver that oneAPI.jl ships (NEO 26.x).

Ready:
- The GPU code is written for any KernelAbstractions backend; `ext/DolooneAPIExt.jl` loads,
  with nothing specific to oneAPI.
- Tests: `test/gpu/gpu_tests.jl`, shared by CUDA (`test/cuda`), oneAPI (`test/oneapi`) and
  JLArrays (in the main test suite). On devices without Float64, Float32 only.
- Benchmarks: `--backend=oneapi`, in the `benchmark/oneapi` environment (see benchmark/README.md).

Checked without an Intel GPU:
- All GPU code paths, solvers included, pass with JLArrays (generic GPU array: no scalar
  indexing, kernels through KernelAbstractions).
- Float32 kernels contain no Float64 instruction (checked in the PTX of every kernel on CUDA).
  The cell index of the splines (`div` of Float32, computed in Float64 by Base) was the only one.

On the laptop:

```bash
git checkout integration-gpu
julia --project=test/oneapi -e 'using Pkg; Pkg.instantiate(); using oneAPI; oneAPI.versioninfo()'
# resources, and what integration_B needs (memory, Float64 support, kernel arguments)
julia --project=benchmark/oneapi -e 'using Pkg; Pkg.instantiate(); using oneAPI, Dolo;
    model = include(joinpath(pkgdir(Dolo), "examples/ymodels/integration_B.jl"));
    display(Dolo.resources(model))'
julia --project=test/oneapi test/oneapi/runtests.jl
julia -t auto --project=benchmark/oneapi benchmark/integration_B.jl --backend=oneapi --sizes=10,20 --engines=threads,gpu
julia -t auto --project=benchmark/oneapi benchmark/integration_B_solve.jl --backend=oneapi --sizes=20
```

To watch for:
- Kernel arguments: the residual kernels receive `(dmodel, x, φ)` by value, 2.0 KB in Float32
  for integration_B (3.1 KB in Float64), more than CUDA's usual kernels. If Intel's limit is
  lower, launches fail: then pass the large constant parts (Markov transition matrix, grids
  of the 9 splines) in device arrays.
- GMRES relies on oneMKL for `dot`, `norm` and `mul!` on views of `oneArray`.
- Memory is shared with the system on integrated GPUs. Workspace in Float32: about
  4,900 bytes per point with GMRES (2,600 for `L`, 1,600 for the Krylov basis), so
  350 MB at 20³ and 5.6 GB at 50³.
- Tolerances in Float32: `tol_ε=1e-4`, `linsolve_tol=1e-5` (the benchmark scripts set them).

## GPU: integration_B on the 50³ grid in Float64, with Neumann

Try solving integration_B at 50³ (1,125,000 points) in Float64 on the GPU with
`linsolve=:neumann`, with the integrated graphics card driving the display so that the
whole discrete GPU is available.

- With the display on the RTX 2080, 7,470 MB of its 8,162 MB are free. The Newton
  workspace with Neumann needs about 7,330 MB, so it allocates (190 MB left), but the
  spline fit then needs a temporary buffer of 191 MB.
- With GMRES it needs 11.2 GB (the Krylov basis alone is 3.7 GB): out of reach either way
  without reducing the memory used.
- Measure with `CUDA.memory_info()` before allocating; components and sizes per grid
  point were measured on 2026-10-19 (`L`: 5,220 bytes/point, of which `F_x` 4,608).

## Other follow-ups from the GPU work

- Memory: store `F_x` in Float32 in Float64 solves (−2.6 GB at 50³), drop `S.val` from
  `BlockVector` (unused when applying `L`, −0.4 GB), Krylov basis in Float32.
- GMRES: Gram-Schmidt with one reduction per iteration (fewer GPU-CPU syncs and passes
  over the basis). Its own work is about 18% of Newton at 20³.
- Applying `L` reaches 70–90 GB/s out of 448 GB/s: look at the reads of the spline
  coefficients.
- Default tolerances by precision: `tol_ε = 1e-8` is out of reach in Float32 (each time
  iteration step then runs its maximum number of Newton steps).
- Newton from the initial guess fails on integration_B (negative number to a fractional
  power); the warm-up time iteration steps cost about as much as Newton.
- Accelerated time iteration gives NaN on `consumption_savings_ar1` (same with the YAML
  version of the model).
