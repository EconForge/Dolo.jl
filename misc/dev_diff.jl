using Dolo
using Dolo: transition

root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/consumption_savings.jl")

dm = Dolo.discretize(model)
wk = Dolo.time_iteration_workspace(dm)


x0 = wk.x0
φ = wk.φ
r0 = Dolo.F(dm, x0, φ)


s = dm.grid[1 + im]
x = x0[1]


# both work
# the second is two times slower and makes a few allocations

using ForwardDiff
ForwardDiff.jacobian(u->Dolo.F(dm, s, u, φ), x)


using Enzyme
Enzyme.jacobian(Forward, u->Dolo.F(dm, s, u, φ), x)


using StaticArrays
function diffit(dm,s,x,φ)
    return sum(Enzyme.jacobian(Forward, u->Dolo.F(dm, s, u, φ), x))
end


using BenchmarkTools
@time diffit(dm, s,x, φ);
@time Dolo.dF_1(dm,s,x,φ);

@benchmark Dolo.dF_1(dm,s,x,φ)
@benchmark diffit(dm, s,x, φ)

@code_llvm diffit(dm, s,x, φ)



# io = IOBuffer()
# code_llvm(io,Dolo.dF_1,[typeof(e) for e in (typeof,s,x,φ)])
# c_1 = String(take!(io))

# io = IOBuffer()
# code_llvm(io,diffit,[typeof(e) for e in (typeof,s,x,φ)])
# c_2 = String(take!(io))


# reverse with Enzyme doesn't work
jacobian(Reverse, u->Dolo.F(dm, s, u, φ), x, Val(2))




using StaticArrays
u0 = SVector(0.2, 0.3)
fun(u::SVector{2,Float64}) = u.^-2
Enzyme.jacobian(Forward, fun,  u0)


Enzyme.jacobian(Reverse, fun,  u0, Val(2))


using Enzyme
Enzyme.onehot(u0)
typeof(Enzyme.onehot(u0))



BatchDuplicated(m::SVector{2,Float64}, t::Tuple{MVector{2,Float64},MVector{2,Float64}}) = 
    BatchDuplicated(MVector(m...), t)



w = SVector(0.2, 0.3)
hun(u::SVector{2,Float64}) = u.^-2
hun(u::MVector{2,Float64}) = gun(SVector(u))
mw = MVector(w...)


res = Enzyme.jacobian(Forward, hun,  mw)



# ERROR: MethodError: no method matching BatchDuplicated(::SVector{2, Float64}, ::Tuple{MVector{2, Float64}, MVector{2, Float64}})

# BatchDuplicated(SVector(0.2, 0.4), (MVector(0.3, 2.0), MVector(0.3, 9.5)))



Enzyme.jacobian(Forward, u->gun(SVector(u[1],u[2])),  [w...])



fun(u::SVector{2,Float64}) = u.^-2
u0 = SVector(1.0, 2.0)

# naive
Enzyme.jacobian(
    Forward,
    fun,
    u0
)
# ERROR: MethodError: no method matching BatchDuplicated(::SVector{2, Float64}, ::Tuple{MVector{2, Float64}, MVector{2, Float64}})

Enzyme.jacobian(
    Forward,
    u->fun(SVector(u)),
    MVector(u0)
)

function doit(u0)
    mat = Enzyme.jacobian(
        Forward,
        u->fun(SVector(u)),
        MVector(u0)
    )
    return sum(mat)
end

@time doit(w)