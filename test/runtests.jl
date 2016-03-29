using Dolo
using Compat
using DataStructures

# write your own tests here
if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

include("parser.jl")
include("util.jl")
include("model_types.jl")
