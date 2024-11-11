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

    import LinearAlgebra: cross, norm, ×

    include("splines/splines.jl")
    using .splines: interp
    

    ⟂(a,b) = min(a,b)
    function ⫫(u,v)
        sq = sqrt(u^2+v^2)
        p =   (v<Inf ? (u+v-sq)/2 : u)
        return p
    end


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

