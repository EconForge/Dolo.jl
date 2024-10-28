using Dolo

root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/rbc_mc.jl")

dm = Dolo.discretize(model)

wksp = Dolo.time_iteration_workspace(dm)

(;x0, φ) = wksp

s = Dolo.QP((1,1),dm.grid[1])
x = x0[1]

f(u) = Dolo.F(dm, s, u, φ)
g(u) = Dolo.F(dm, s, SVector{2,Float64}(u[1],u[2]), φ)


f(x)

function df(u)

    m = MVector(u...)
    fun = u->Dolo.F(dm, s, SVector(u...), φ)
    
    dr = Enzyme.jacobian(Forward, fun, m)
    return dr

end

@time df(x)

@time Dolo.dF_1(dm, s, x, φ)



f(MVector(x...))

function funzyme(u)



# import Zygote
# Zygote.jacobian(f, x)
# Zygote.jacobian(g, [x...])

# import Enzyme
# methods(Enzyme.jacobian)


using Enzyme  # dev enzyme
import Enzyme: BatchDuplicated
using StaticArrays

# function BatchDuplicated(m::SVector{2, Float64}, res::Tuple{MVector{2, Float64}, MVector{2, Float64}})
#     # a = MVector{2,Float64}(m[1],m[2])
#     # b = MVector{2,Float64}(m[1],m[2])
#     # return (a,b)
#     res[1][1] = m[1]
#     res[1][2] = m[2]
#     res[2][1] = m[1]
#     res[2][2] = m[2]
# end

Enzyme.jacobian(Forward,f, x)

Enzyme.jacobian(Forward,g, [x...])

import ForwardDiff

@time ForwardDiff.jacobian(f,x)
# @time Enzyme.jacobian(Reverse, g, [x...])



using StaticArrays
using ForwardDiff


x = SVector(0.1, 0.2)
function test(u)
    return SVector(u[1], u[1]+u[2]^2)
end

x = SVector(0.1, 0.2)
function test2(u)
    return MVector(u[1], u[1]+u[2]^2)
end

test2(x)

e = MVector(x)

function test_mvec(x)

    m = MVector(x...)
    g = Enzyme.jacobian(Forward, test, m)
    
    return sum(g)

end

function test_fd(x)

    g = ForwardDiff.jacobian(test, x)
    
    return sum(g)
end



@time test_mvec(x)

@time test_fd(x)



# function BatchDuplicated(m::SVector{2, Float64}, res::Tuple{MVector{2, Float64}, MVector{2, Float64}})
#     # a = MVector{2,Float64}(m[1],m[2])
#     # b = MVector{2,Float64}(m[1],m[2])
#     # return (a,b)
#     res[1][1] = m[1]
#     res[1][2] = m[2]
#     res[2][1] = m[1]
#     res[2][2] = m[2]
#     nothing
# end


function test_array(u)
    [test(SVector(u[1], u[2]))...]
end

Enzyme.jacobian(Forward, test, x)

Enzyme.jacobian(Forward, test_array, [x...])
