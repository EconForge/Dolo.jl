# Julia-side push

The goal is to have all functionalities from the Python side, while greatly
enhancing the structure. In the meantime, basic compatibility is maintained on
the Python side, so that one can still solve the models.

At this stage we limit the available options to the minimum.
The following items are available to work on:

- AR1 decision rule

    Decision rules are defined on processes and/or discretized processes.

- Simulate methods for all processes

    `simulate(process; start=Vector{Float64}, T=40)::DataFrame`

- Simulate a model with a solution

    `simulate(model, dr; start=Vector{Float64}, T=40)::DataFrame`

- Solution algorithms

    - VFI (cleanup)
    - PEA
    - perturbation (requires fix to AR1)

- Distribution object:

    Here is how the type could look like:

        `Distribution(grid_endo, grid_exo, density)`   # sum(density)==1
        `Distribution(grid, density)`

    It can be constructed from observations by assigning each observation to the
    nearest gridpoint.

- Ergodic distribution:

    Given a model and a decision rule, one can compute the ergodic distribution
    of the states, either by simulating, or by solving a markov chain describing
    the transitions.

- Accuracy measures

# Cleanup and documentation

Define which are the right options to pass to all methods, in a consistent and
predictible fashion. (i.e. options for outer loop, for inner loop, for endo grid, etc.)
Avoid extreme julianisms like `BSpline(Quadratic(Reflect())), OnCell()` but maintain correct performances.

Document all user facing methods. Write a description for each algorithm (language independent).

Move all examples in a separate repository (yaml and notebooks).

# Improvements

## Language

Implement required changes in the yaml specification.

## API and internals

- add all kind of derivative to decision rules (w.r.t. states, coefficients, ...)
- use custom sparse matrices objects instead of 3d arrays
- use NLSolve.jl instead of basic newton (requires fix to NLSolve)
- allow for many types of interpolation (with Interpolations.jl or BasisMatrices.jl)

## Performance

Right now, the Julia methods are faster than Python for smaller models, and lag
behind for big ones. This is *not* a problem until there is an initial implementation
for everything. Then there are many possible avenues:

- improve vectorized evaluations of the model
- perform inplace computations where it matters
- implement a decision rule which precomputes the bases for tomorrows exogenous values (for AR1)
- use derivative information in some of the methods (e.g. one step value evaluation)
