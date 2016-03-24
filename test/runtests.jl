using Dolo
using Base.Test

# write your own tests here
if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

include("parser.jl")
