# Algorithms

```@meta
CurrentModule = Dolo
```
## Time iteration algorithm:

```@docs
time_iteration
```

Global
- verbose::Bool=true
- details::Bool=true
- interpolation (multilinear, cubic, chebychev, ...)
    - ...
- discretization:
    - `type::Symbol=rouwenhorst`
    - options: (rouwenhorst)
        - `N::Int=3`
    - options: (tauchen)
        - `N::Int=3`
        - `mu::Int=3`, number of standard deviations
    - options: (gdp)
        - `mu::Float=3`, number of standard deviations
        - `N::5`, number of points today
        - `N_int::5` number of integration nodes in each dimension
- outer loop control:
    - `maxit::Int=100`
    - `tol_ϵ::Float=1e-8`
    - `tol_η::Float=(1-λ)*tol_ϵ`
    - how to estimate lambda (def: 0.99)
- complementarities:
    - `complementarities::Bool=true`
- inner loop control:
    - `solver-type= safeguarded-newton`
    - solver-options: (for newton)
        - `verbose::Bool=false`
        - `maxit::Int=10`
        - `tol::Float=1e-6`
        - `eps::Float=1e-8` (to evaluate numerical jacobian)
        - `n_bsteps::Int=5`
        - `lam_bsteps::Float=0.5`
- whether to use expectation function
    - `use_expectations::Bool=false` (true requires functions `fh`)
    - `use_direct_response::Bool=false` (uses explicit formula to to solve arbitage equations, requires `dh`)

Full signature (for now):

```
time_iteration(model, discretized_process, endogenous_grid, initial_guess;
    maxit, tol_ϵ, tol_η, λ,
    complementarities,
    grid=Dict(...), # overrides get_grid() values
    solver=Dict(...),  # options for solver, unpacked an passed to solver
    domain=Dict(...), # overrides get_domain() values
    use_expectations, use_direct_response
)
```

## Time iteration algorithm (direct):

```@docs
time_iteration_direct
```


## Value iteration

```@docs
value_iteration
```


Options:

- verbose::Bool=true
- details::Bool=true
- interpolation (multilinear, cubic, chebychev, ...)
    - ...
- discretization:
    - `method::Symbol=rouwenhorst`
    - `N::Int=3`
- outer loop control:
    - `maxit::Int=100`
    - `tol_ϵ::Float=1e-8`
    - `tol_η::Float=(1-λ)*tol_ϵ`
    - how to estimate lambda (def: 0.99)
- complementarities:
    - `complementarities::Bool=true`
- inner loop control:
    - `solver-type= safeguarded-newton`
    - solver-options:
        - `verbose::Bool=false`
        - `maxit::Int=10`
        - `tol::Float=1e-6`
        - `eps::Float=1e-8` (to evaluate numerical jacobian)
        - `n_bsteps::Int=5`
        - `lam_bsteps::Float=0.5`
- whether to use expectation function
    - `use_expectations::Bool=false` (true requires functions `fh`)
    - `use_direct_response::Bool=false` (uses explicit formula to to solve arbitage equations, requires `dh`)
