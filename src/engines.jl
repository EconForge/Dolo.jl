# Execution engines
#
# Every vectorized operation on a grid is written as a *pointwise kernel*
#
#     f(n, args...) -> value at grid point n
#
# which is applied to all grid points by `pointwise!(f, out, engine, args...)`.
# Kernels are the canonical (non-allocating) implementation of an operation,
# engines only decide how the loop over grid points is run.

import KernelAbstractions
import KernelAbstractions: @kernel, @index

abstract type Engine end

"Plain loop on one CPU thread."
struct Serial <: Engine end

"Loop split among Julia threads (`Threads.@threads`)."
struct Threaded <: Engine end

"KernelAbstractions.jl kernel, on any backend (GPU, or `KernelAbstractions.CPU()`)."
struct KA{B} <: Engine
    backend::B
end

"""
    pointwise!(f, out, engine, args...)

Set `out[n] = f(n, args...)` for every (linear) index `n` of `out`.
(`F` and `N` force specialization on `f` and on the arguments.)
"""
function pointwise! end

function pointwise!(f::F, out, ::Serial, args::Vararg{Any,N}) where {F,N}
    for n in 1:length(out)
        @inbounds out[n] = f(n, args...)
    end
    return out
end

function pointwise!(f::F, out, ::Threaded, args::Vararg{Any,N}) where {F,N}
    Threads.@threads for n in 1:length(out)
        @inbounds out[n] = f(n, args...)
    end
    return out
end

@kernel function pointwise_kernel!(out, f, args)
    n = @index(Global, Linear)
    @inbounds out[n] = f(n, args...)
end

function pointwise!(f::F, out, engine::KA, args::Vararg{Any,N}) where {F,N}
    pointwise_kernel!(engine.backend)(out, f, args; ndrange=length(out))
    KernelAbstractions.synchronize(engine.backend)
    return out
end


"""
    select_engine(engine, x)

Engine from a user option: an `Engine`, or one of `:serial`, `:threads`, `:ka`
(KernelAbstractions on the backend where the data `x` lives, e.g. a GPU).
"""
select_engine(engine::Engine, x) = engine

function select_engine(engine::Union{Symbol,Nothing}, x)
    if engine in (:serial, :none, nothing)
        return Serial()
    elseif engine in (:threads, :cpus, :cpu)
        return Threaded()
    elseif engine in (:ka, :gpu)
        return KA(KernelAbstractions.get_backend(x.data))
    else
        error("Unknown engine $engine. Use :serial, :threads or :ka.")
    end
end
