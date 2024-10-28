using Dolo.splines
using StaticArrays

ranges = (
    (0.0, 1.0, 10),
    (0.0, 1.0, 10),
    (0.0, 1.0, 10)

)
a = tuple((r[1] for r in ranges)...)
b = tuple((r[2] for r in ranges)...)
n = tuple((r[3] for r in ranges)...)
c = [sin(i+j+k) for i=1:12, j=1:12, k=1:12]


s0 = SVector(-.5,-0.5-0.5)
s1 = SVector(0.5,0.5,0.5)
s2 = SVector(1.5,1.5,1.5)


import Dolo
using BenchmarkTools

@time res1 = Dolo.splines.eval_UC_spline(a,b,n,c,s1)

@time res = Dolo.splines.eval_spline(ranges,c,s1)

@code_llvm Dolo.splines.eval_spline(ranges,c,s1)




A = Dolo.splines.mextract(c, (1,1,1),Val(4))
Φ = tuple( (Dolo.splines.cPhi(e) for e in (0.5,0.5,0.5))... )


@code_llvm Dolo.splines.reduce_tensors(A,Φ)

@benchmark res1 = Dolo.splines.eval_UC_spline(a,b,n,c,s2)
@benchmark res = Dolo.splines.eval_spline(ranges,c,s2)

@time res2 = Dolo.splines.eval_spline(ranges,c,s)


Dolo.splines.create_Phi(2,"natural",false)


function mextract(c::AbstractArray{T,d}, inds::NTuple{d,<:Int}, ::Val{n}) where T where d where n
    ranges = tuple((i:i+n-1 for i in inds)...)
    return SArray{NTuple{d,n},T,d,d^n}(view(c, ranges...))
end