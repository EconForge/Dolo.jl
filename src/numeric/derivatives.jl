"""
```julia
eval_jacobian(m::ANM, f!::Function, n::Int, before::Tuple, to_diff::Tuple,
              after::Tuple, i::Int)
```

Evaluate the Jacobian for the `f!` function of model `m`. `n` is the number
of equations in `f!`. `before` is a tuple of arguments to `f!` that come before
the argument to be differentiated. `to_diff` is a tuple of the arguments that
might be differentiated, `after` is a tuple of arguments that come after the one
to be differentiated, and `i` specifies which element of `to_diff` should be
differentiated. Thus, `f!` is called like this:

```julia
f!(out, m, before..., to_diff[1:i-1]..., to_diff[i] ,to_diff[i+1:end]...,
   after...)
```
"""
function eval_jacobian(m::ANM, f!::Function, n::Int, before::Tuple,
                       to_diff::Tuple, after::Tuple, i::Int)
    _f!(out, _) = f!(m, out, before..., to_diff[1:i-1]..., _,
                     to_diff[i+1:end]..., after...)
    j! = ForwardDiff.jacobian(_f!; mutates=true, output_length=n)
    out = Array(Float64, n, length(to_diff[i]))
    j!(out, to_diff[i])
    out
end

"""
```julia
eval_jacobians(m::ANM, f!::Function, n::Int, before::Tuple, to_diff::Tuple,
               after...)
```

Compute certain Jacobians of the function `f!` for model `m`. `n` is the number
of equations in `f!`. `before` is a tuple of arguments before those to be
differentiated, `to_diff` is a tuple of arguments to differentiated, `after`
is all arguments that come in `f!` after the arguments for the derivatives.

The output is a vector containing `length(to_diff)` `Matrix{Float64}` elements.
"""
function eval_jacobians(m::ANM, f!::Function, n::Int, before::Tuple,
                        to_diff::Tuple, after...)
    out = Matrix{Float64}[eval_jacobian(m, f!, n, before, to_diff, after, i)
                          for i in 1:length(to_diff)]
end

"""
```julia
eval_jacobians(m::ANM, f!::Function, n::Int, to_diff::Tuple, after...)
```

Version where no arguments come `before` `to_diff` when calling `f!`. The
full version of `eval_jacboians` is called with `before` set to an empty
tuple
"""
eval_jacobians(m::ANM, f!::Function, n::Int, to_diff::Tuple, after...) =
    eval_jacobians(m, f!, n, tuple(), to_diff, after...)
