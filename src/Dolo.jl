__precompile__(false)
module Dolo

    # include("./stupid.jl")

    import Base: eltype
    # using LabelledArrays
    using StaticArrays
    using ForwardDiff
    using Printf
    using Format 
    using Crayons.Box
    using Term

    import LinearAlgebra: cross, norm, ×

    include("splines/splines.jl")
    using .splines: interp
    

    ⟂(a,b) = min(a,b)

    function ⫫(u0,v0)
        BIG = 100000
        u = min(max(u0,-BIG),BIG)
        v = min(max(v0,-BIG),BIG)
        sq = sqrt(u^2+v^2)
        p =   (u+v-sq)/2
        return p
    end

    # function ⫫(u,v)
    #     sq = sqrt(u^2+v^2)
    #     p =   (v<Inf ? (u+v-sq)/2 : u)
    #     return p
    # end

    using ChainRulesCore
    using ForwardDiffChainRules

 # define your frule for function f1 as usual
function ChainRulesCore.frule((_, Δu, Δv), ::typeof(⫫), u::Real, v::Real)
    BIG = 100000
    u = min(max(u0,-BIG),BIG)
    v = min(max(v0,-BIG),BIG)
    sq = sqrt(u^2+v^2)
    Omega = (u+v-sq)/2
    if u==v==0.0
        Omega, (Δu + Δv)*0
    # elseif !(u<Inf)
    #     return v, Δv
    # elseif !(v<Inf)
    #     return u, Δu
    else
        Δ1 = (0.5 - u/sq)*Δu
        Δ2 = (0.5 - v/sq)*Δv
        return Omega,  Δ1 + Δ2
    end
 end


# @ForwardDiff_frule ⫫(u::ForwardDiff.Dual, v::ForwardDiff.Dual)
# @ForwardDiff_frule ⫫(u::ForwardDiff.Dual, v)
# @ForwardDiff_frule ⫫(u, v::ForwardDiff.Dual)


# ⫫(u,v) = fun(u,v)

    #   #, SDiagonal(J_u), SDiagonal(J_v)
    #     # J_u = (v<Inf ? (1.0 - u[i]./sq[i])/2 : 1) for i=1:d )
    #     # J_v = (v<Inf ? (1.0 - v[i]./sq[i])/2 : 0) for i=1:d )
    
    #     return p  #, SDiagonal(J_u), SDiagonal(J_v)
    # end

    import Base: getindex

    # import LabelledArrays: merge

    # TODO
    converged(sol::NamedTuple) = (sol.message == "Convergence")


    # # type piracy
    # function merge(a::SLArray, b::NamedTuple)
    #     @assert issubset(keys(b), keys(a))
    #     SLVector( (merge(NamedTuple(a), b)) )
    # end


    struct QP{A,B}
        loc::A
        val::B
    end

    import Base: show
    
    show(io::IO, x::QP) = print(io,"QP(;loc=$(x.loc), val=$(x.val))")
    
    macro I(v)
        :(($v).index)
    end

    macro V(v)
        :(($v)[2])
    end

    

    include("misc.jl")
    include("space.jl")
    include("grids.jl")   
    include("processes.jl")
    include("garray.jl")
    include("model.jl")
    include("funs.jl")
    include("conversion.jl")

    include("dev_L2.jl")
    include("algos/simul.jl")
    include("algos/results.jl")
    include("algos/time_iteration.jl")
    include("algos/time_iteration_accelerated.jl")
    include("algos/vfi.jl")


    include("adapt.jl")
    include("utils.jl")

    # function yaml_import(filename)
    #     DoModel.DoloModel(filename)
    # end

    export time_iteration, tabulate, discount_factor

    module Build
        using Dolo: transition, arbitrage, bounds, reward
        using Dolo: GSpace, CSpace
        using Dolo: SGrid, CGrid
        using Dolo: CartesianSpace, GridSpace
        using Dolo: MvNormal, MarkovChain
        using Dolo: YModel
        export SGrid, CGrid, CSpace, GSpace, CartesianSpace, GridSpace, transition, arbitrage, MvNormal, MarkovChain
        export YModel, bounds, reward
        export rouwenhorst
        export SMatrix, SVector
    end


end # module

