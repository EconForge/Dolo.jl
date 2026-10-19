# TODO

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
