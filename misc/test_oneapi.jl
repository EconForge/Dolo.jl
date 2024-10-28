using oneAPI

using KernelAbstractions
using KernelAbstractions: synchronize


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
# using CUDA
import Adapt: adapt_structure, adapt
import Adapt

# using CUDAKernels
using KernelAbstractions
using KernelAbstractions: @kernel, @index, @Const
# using CUDA


adapt_structure(to, s::Dolo.SplineInterpolator{G,C,k}) where G where C where k = let
  θθ = adapt(to, s.θ)
  CC = typeof(θθ)
  Dolo.SplineInterpolator{G,CC,k}(s.grid, θθ)
end


T = Float32
TN = Int64

M = 100
grid = (
    (convert(T,0.0),convert(T,1.0),convert(TN,M)),
    (convert(T,0.0),convert(T,1.0),convert(TN,M)),
    (convert(T,0.0),convert(T,1.0),convert(TN,M))
)


C = rand(T, M+2,M+2,M+2)

spl = Dolo.SplineInterpolator{typeof(grid), typeof(C), 3}(grid, C)

gspl = adapt(oneArray, spl)

cg = Dolo.CGrid(grid)

# gcg = adapt(cg)
res = [spl(e) for e in cg];
res32 = convert(Vector{Float32}, res*0);
gres32 = adapt(oneArray,res32);


using StaticArrays

@kernel function mykernel(grid,r,sp)

  I = @index(Global)
  p = grid[I]
  r[I] = sp(p)

end 

function time_kernel(grid, res, sp; N=64, K=1000)
  backend = get_backend(res)
  for k=1:K
    mykernel(backend, N)(grid, res, sp, ndrange=size(grid))
  end
  synchronize(backend)
  # rr = adapt(Array,gres32)
end
  # rr = adapt(Array,gres32)

@time time_kernel(cg, res32, spl; N=16)
@time time_kernel(cg, res32, spl; N=32)
@time time_kernel(cg, res32, spl; N=64)
@time time_kernel(cg, res32, spl; N=128)
@time time_kernel(cg, res32, spl; N=256)


@time time_kernel(cg, gres32, gspl; N=16)
@time time_kernel(cg, gres32, gspl; N=32)
@time time_kernel(cg, gres32, gspl; N=64)
@time time_kernel(cg, gres32, gspl; N=128)
@time time_kernel(cg, gres32, gspl; N=256)


@time rres = [spl(e) for e in cg];

rr = adapt(Array,gres32)
maximum(abs.(rr - rres))
# import Base: floor

# Base.floor(a::Float32) = a <1.0 ? 0 : 10

# trunc(::Type{Signed}, x::IEEEFloat) = trunc(Int,x)
# trunc(::Type{Unsigned}, x::IEEEFloat) = trunc(UInt,x)
# trunc(::Type{Integer}, x::IEEEFloat) = trunc(Int,x)

# # fallbacks
# floor(::Type{T}, x::AbstractFloat) where {T<:Integer} = trunc(T,round(x, RoundDown))
# ceil(::Type{T}, x::AbstractFloat) where {T<:Integer} = trunc(T,round(x, RoundUp))
# round(::Type{T}, x::AbstractFloat) where {T<:Integer} = trunc(T,round(x, RoundNearest))



synchronize(ev)


backcpu = get_backend(res32)
ev = mul2_kernel(backcpu, 64)(cg,res32,spl, ndrange=size(cg))



import Dolo: splines

splines.create_local_parameters(1; Tf=Float32)

splines.create_Phi(2, "natural", false; Tf=Float32)

splines.create_function(1,"natural"; vectorize=false,Tf=Float32)


fun = include("p.jl")

CC = rand(Float32, 10, 10, 10)
using StaticArrays
p = SVector(0.5f0, 0.5f0, 0.5f0)

eval_UC_spl( (0.0f0, 0.0f0, 0.0f0), (1.0f0, 1.0f0, 1.0f0), [8,8,8], CC, p)



using oneAPI
using KernelAbstractions
using Adapt
const KAA = KernelAbstractions

@kernel function round_kernel(A,B)
  I = @index(Global)
  v = trunc(B[I])
  A[I] = (v::Int32)
end


A = zeros(Int8, 1024, 1024)
B = rand(Float32, 1024, 1024)*10
A_gpu = adapt(oneArray, A)
B_gpu = adapt(oneArray, B)
dev = get_backend(A_gpu)
round_kernel(dev, 64)(A_gpu, B_gpu, ndrange=size(A_gpu))
KAA.synchronize(dev)
