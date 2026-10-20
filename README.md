# Dolo

![Build Status](https://github.com/EconForge/Dolo.jl/workflows/CI/badge.svg)


Documentation: [http://www.econforge.org/Dolo.jl/latest/](http://www.econforge.org/Dolo.jl/latest/)

Dolo is a tool to describe and solve economic models. It provides a simple classification scheme to describe many types of models, represents them as efficient Julia objects, and provides many reference solution algorithms to find the solution of these models under rational expectations.

Models are written directly in Julia (see [examples/ymodels](examples/ymodels)). To write them as simple YAML text files instead, use [DoloYAML.jl](https://github.com/EconForge/DoloYAML.jl), which imports them as models fully compatible with Dolo.jl.

Dolo understand several types of nonlinear models with occasionally binding constraints (with or without exogenous discrete shocks), as well as local pertubations models.

## Diagnostics

`Dolo.diagnose()` reports the computing resources available to Dolo: threads, memory, and GPUs
(when their package is loaded: `using CUDA` or `using oneAPI`). Given a model, or a discretized
model, it also reports what solving it needs, and whether it fits:

```julia
using Dolo, CUDA
model = include("examples/ymodels/integration_B.jl")
Dolo.diagnose(model; discretization=Dict(:endo => [50, 50, 50]))
```

- size of the grid, memory of each solver (time iteration, Newton with Neumann or GMRES) in the
  precision of the model and in Float32, compared with the memory of the CPU and of each GPU;
- size of the arguments of the GPU kernels;
- checks: the discretized model can be passed to GPU kernels, the residuals at the initial guess
  are finite, a Float32 copy of the model exists and computes in Float32;
- recommendations (number of threads, precision on the GPU, linear solver).

More docs coming soon...
