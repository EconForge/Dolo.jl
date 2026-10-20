module Dolo

    import Base: eltype
    using StaticArrays
    using ForwardDiff
    using Printf
    using Format 
    using Crayons.Box
    using Term

    import LinearAlgebra: cross, norm, ×

    include("splines/splines.jl")
    using .splines: interp
    # the cell of a dual number is the cell of its value
    splines.cell_index(u::ForwardDiff.Dual, δ) = splines.cell_index(ForwardDiff.value(u), δ)
    

    ⟂(a,b) = min(a,b)

    # Fischer-Burmeister complementarity (arguments are clamped so that infinite bounds are allowed)
    function ⫫(u0,v0)
        BIG = 100000
        u = min(max(u0,-BIG),BIG)
        v = min(max(v0,-BIG),BIG)
        sq = sqrt(u^2+v^2)
        p =   (u+v-sq)/2
        return p
    end

    using ChainRulesCore
    using ForwardDiffChainRules

    # derivative of ⫫, well defined at u==v==0
    # not registered with ForwardDiff yet (see @ForwardDiff_frule)
    function ChainRulesCore.frule((_, Δu0, Δv0), ::typeof(⫫), u0::Real, v0::Real)
        BIG = 100000
        u = min(max(u0,-BIG),BIG)
        v = min(max(v0,-BIG),BIG)
        Δu = (-BIG < u0 < BIG) ? Δu0 : zero(Δu0)
        Δv = (-BIG < v0 < BIG) ? Δv0 : zero(Δv0)
        sq = sqrt(u^2+v^2)
        Omega = (u+v-sq)/2
        if u==v==0
            return Omega, (Δu + Δv)*0
        else
            Δ1 = (1 - u/sq)/2*Δu
            Δ2 = (1 - v/sq)/2*Δv
            return Omega,  Δ1 + Δ2
        end
    end

    import Base: getindex

    # TODO
    converged(sol::NamedTuple) = (sol.message == "Convergence")


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

    include("engines.jl")
    include("algos/linsolve.jl")
    include("linear_operator.jl")
    include("algos/residuals.jl")
    include("algos/simul.jl")
    include("algos/results.jl")
    include("algos/time_iteration.jl")
    include("algos/vfi.jl")


    include("adapt.jl")
    include("utils.jl")
    include("resources.jl")

    export time_iteration, newton, tabulate, discount_factor

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

