# test/runtests.jl
using Dolo
using Test

@testset "General" begin
    include("general.jl")
end