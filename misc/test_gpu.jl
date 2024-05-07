using Dolo
using StaticArrays
using BenchmarkTools
import Dolo.splines: CubicInterpolator
using LoopVectorization
# using CUDAKernels
using KernelAbstractions
using KernelAbstractions: synchronize
using KernelAbstractions: @kernel, @index, @Const
using Adapt: adapt
using CUDA
import Adapt: adapt_structure, adapt
import Adapt

# using CUDAKernels
using KernelAbstractions
using KernelAbstractions: @kernel, @index, @Const
using CUDA


grid = (
    (0.0,1.0,20),
    (0.0,1.0,20),
    (0.0,1.0,20)
)


# C = rand(SVector{2, Float64}, 22,22,22)
C = rand(22,22,22)

spl = Dolo.SplineInterpolator{typeof(grid), typeof(C), 3}(grid, C)

# function test(spl; N=100000)

#     out = zeros(eltype(spl.θ), N)
#     for n=1:N
#         t = sin(n)/2.0+1.0
#         x = SVector(t,1-t,t^2)
#         r = spl(x)
#         out[n] = r
#     end

#     return sum(out)

# end


# @benchmark  test(spl)


# function test2(spl; N=100000)

#     out = zeros(N)
#     @tturbo for n=1:N
#         t = sin(n)/2.0+1.0
#         x = SVector(t,1-t,t^2)
#         r = spl(x)
#         out[n] = r
#         # out[n] = r
#     end

#     return sum(out)

# end

# @benchmark test2(spl)

# # using CUDAKernels

# @kernel function ker_(out,(spl))
#     n = @index(Global)
#     t = sin(n)/2.0+1.0
#     x = SVector(t,1-t,t^2)
#     r = spl(x)
#     out[n] = r
# end

# function test3(spl; N=100000)

#     out = zeros(N)

#     # backend = get_backend(out)

#     kernel! = ker_(CPU())

#     done = kernel!(out, spl, ndrange=N)
#     wait(done)

#     return sum(out)

# end

# @benchmark test3(spl)


# function test4(spl; N=100000)

#     out = CUDA.zeros(eltype(spl.θ), N)

#     kernel! = ker_(CUDADevice())

#     done = kernel!(out, spl, ndrange=N)
#     wait(done)

#     return sum(out)

# end


adapt_structure(to, s::Dolo.SplineInterpolator{G,C,k}) where G where C where k = let
    θθ = adapt(to, s.θ)
    CC = typeof(θθ)
    Dolo.SplineInterpolator{G,CC,k}(s.grid, θθ)
end


gspl = adapt(CuArray, spl)

# test4(gspl)

# @benchmark test4(gspl)


function orphan(model::DoloYAML.Model{ID, Dom}) where ID where Dom
    P = typeof(model.parameters)
    DoloYAML.DummyModel{ID, Dom, P}(model.parameters)
end

function orphan(model::Dolo.YModel)
    Dolo.YModel(
        Dolo.name(model),
        model.states,
        model.controls,
        model.exogenous,
        model.calibration,
        orphan(model.source)
    )
end


function orphan(model::Dolo.DYModel)
    Dolo.DYModel(
        orphan(model.model),
        model.grid,
        model.dproc
    )
end



root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/rbc_mc.jl")


dmodel = Dolo.discretize(model)



# model2 = orphan(model)
# dmodel2 = orphan(dmodel)


x0 = Dolo.GVector(dmodel.grid, [Dolo.calibrated(model, :controls) for i=1:length(dmodel.grid)])

φ = Dolo.DFun(model, x0)

Dolo.F(dmodel, x0, φ)

Dolo.calibrated(model, :controls)

import Adapt: adapt_storage

function adapt_structure(to, v::Dolo.GVector{G,V}) where G where V
    data = adapt(to, v.data)
    Dolo.GVector{G, typeof(data)}(
        v.grid,
        data
    )
end
function adapt_structure(to, f::Dolo.DFun{Dom, Gar, Itp, vars}) where Dom where Gar where Itp where vars
    itp = adapt(to, f.itp)
    values = adapt(to, f.values)
    Dolo.DFun{Dom, typeof(values), typeof(itp), vars}(
        f.domain,
        values,
        itp
    )
end

@kernel function ker_F!(out, dmodel, x0, φ)

    n = @index(Global, Linear)

    i,j = Dolo.from_linear(dmodel.grid, n)
    # I = @index(Global, Cartesian)

    s0 = Dolo.QP( (i,j), dmodel.grid[n])
    x = x0[n]

    r = Dolo.F(dmodel, s0, x, φ)

    out[n] = r

    # @synchronize
end

# tt = typeof(model)


# struct MyType{A,B}
#     a::A
#     b::B
# end

# mm = MyType((1,2,3),(:a, :b))


# source = :(MyType{Tuple{Int64, Int64, Int64}, Tuple{Symbol, Symbol}})
# dest = :(MyType{Tuple{Int64, Int64, Int64}})


# t_source = eval(source)
# t_dest = eval(dest)


function F_cpu(dmodel, x0, φ, K=1000)

    backend = CPU()

    r0 = x0*0                     # TODO

    # kernel! = ker_F(CUDADevice())
    kernel! = ker_F!(backend)

    # N = length(dmodel.grid)
    N1 = length(dmodel.grid.g1)
    N2 = length(dmodel.grid.g2)
    for k=1:K
        kernel!(r0, dmodel, x0, φ; ndrange=(N1,N2))
    end

    # @synchronize
    r0 = adapt(Array, r0)

    return sum(sum(r0.data))
    
end

@time res = F_cpu(dmodel, x0, φ);



function F_gpu(dmodel, x0, φ, K=1000)
    
    backend = CUDABackend()
    r0 = x0*0                     # TODO
    r0_gpu = adapt(CuArray, r0)   # TODO
    x0_gpu = adapt(CuArray, x0)   # TODO
    φ_gpu = adapt(CuArray, φ)

    # kernel! = ker_F(CUDADevice())
    kernel! = ker_F!(backend)

    # N = length(dmodel.grid)
    N1 = length(dmodel.grid.g1)
    N2 = length(dmodel.grid.g2)

    for k=1:K
        kernel!(r0_gpu, dmodel, x0_gpu, φ_gpu, ndrange=(N1,N2))
    end
    r0 = adapt(Array, r0_gpu)
    return sum(sum(r0))

end

@time res1 = F_cpu(dmodel, x0, φ);
@time res2 = F_gpu(dmodel, x0, φ);








function F_fast!(r0, dmodel, x0, φ)

    dm = orphan(dmodel)
    backend = KernelAbstractions.get_backend(x0.data)

    @kernel function ker_FF!(out, dmodel, x0, φ)

        n = @index(Global, Linear)

        i,j = Dolo.from_linear(dmodel.grid, n)
        
        s0 = Dolo.QP( (i,j), dmodel.grid[n])
        x = x0[n]
        r = Dolo.dF_1(dmodel, s0, x, φ)
        out[n] = r
    end

    kernel! = ker_FF!(backend)

    N = length(dm.grid)
    kernel!(r0, dm, x0, φ, ndrange=N)

end



r0 = Dolo.dF_1(dmodel, x0, φ)
r0_gpu = Adapt.adapt(CuArray,r0)
φ_gpu = Adapt.adapt(CuArray,φ)
x0_gpu = Adapt.adapt(CuArray,x0)


function timeit(r0, φ, x0; K=1000)
    for k=1:K
        F_fast!(r0, dmodel, x0, φ)
    end
    sum( sum(r0.data) )
end

@time timeit(r0, φ, x0; K=1)
@time timeit(r0_gpu, φ_gpu, x0_gpu; K=1)

@time sum(r0_gpu.data)


function timeit_FF(r0, φ, x0; K=1000)
    for k=1:K
        FF!(r0, dmodel, x0, φ)
    end
    sum( sum(r0.data) )
end

@time timeit_FF(r0, φ, x0; K=1)

@profview timeit_FF(r0, φ, x0; K=1)
