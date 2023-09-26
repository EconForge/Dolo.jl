# abstract model
abstract type AbstractModel end
abstract type AbstractDModel end

const AModel = AbstractModel
const ADModel = AbstractDModel

macro akin(model)
    typ = typeof(eval(model))
    ntyp = Meta.parse(repr(typ))
    ntyp.args = ntyp.args[1:end-1]
    return ntyp
end

# abstract type AModel{A,B,C,N} end
# abstract type ADModel{A,B,C,D,N} end


# struct Model{A,B,C,N} <: AModel{A,B,C,N}
#     calibration::A
#     endo_domain::B
#     exogenous::C
# end


include("dolo_model.jl")

# struct Model{A,B,C,N} <: AModel
#     calibration::A
#     domain::B
#     transition::C
# end



# struct Model{A,B,C,N} <: AModel{A,B,C,N}
#     calibration::A
#     domain::B
# end

# exo_transition(model, ... )
# endo_transition(model, ...)


# struct Model{A,B,C,N} <: AModel{A,B,C,N}
#     calibration::A
#     domain::B
#     transition::C
# end


# Model(calibration::A, domain::B, exogenous::C; name::Symbol=:anonymous) where A where B where C = Model{A,B,C,name}(calibration, domain, transition)
# name(::Model{A,B,C,N}) where A where B where  C where N = N

# struct DModel{A,B,C,D,N} <: ADModel
#     calibration::A
#     domain::B
#     grid::C
#     transition::D
# end

# DModel(calibration::A, domain::B, grid::C, transition::D; name::Symbol=:anonymous) where A where B where C where D = DModel{A,B,C,D,name}(calibration, domain, grid, transition)

# name(::ADModel{A,B,C,D,N}) where A where B where  C where D where N = N

# function Base.show(io::IO, m::AModel) 
#     println("Model")
#     println("* name = ", name(m))
# end

# function Base.show(io::IO, m::ADModel) 
#     println("Discretized Model")
#     println("* name = ", name(m))
# end

# function recalibrate()
# end

# # TODO
# import Base.show
# Base.show(io::IO, dmodel::ADModel) = print(io, "DModel(#$(hash(typeof(dmodel))))")


# exo_transition(model::ADModel) = model.transition



label_GArray(m, g::GArray) = GArray(g.grid, [LVectorLike(m, e) for e in g.data])

include("model_extensions.jl")

# function arbitrage(model::ADModel, m, s, x, M, S, X, p)
#     m = LVectorLike(model.calibration.m, m)  # this does not keep the original type
#     s = LVectorLike(model.calibration.s, s)
#     x = LVectorLike(model.calibration.x, x)
#     M = LVectorLike(model.calibration.m, M)
#     S = LVectorLike(model.calibration.s, S)
#     X = LVectorLike(model.calibration.x, X)
#     r = arbitrage(model,m,s,x,M,S,X,p)
#     return SVector(r...)
# end

# function arbitrage(model::ADModel, s, x, S, X)
#     p = model.calibration.p
    
#     arbitrage(model, 
#         Dolo.split_states(model, s)...,
#         x,
#         Dolo.split_states(model, S)...,
#         X,
#         p
#     )
# end


# function split_states(model::ADModel, S::Tuple{ind, v}) where ind where v<:SVector
#     split_states(model, S[2])
# end

# function split_states(model::AModel, s_)

#     n_m = ndims(model.exogenous)  # this does not keep the original type
#     n_s = ndims(model.states) - n_m

#     m = SVector((s_[i] for i=1:n_m)...)
#     s = SVector((s_[i] for i=n_m+1:(n_m+n_s))...)

#     return (;m,s)

# end

# TODO: restore initial_guess mechani

# ## default implementation
function initial_guess(model::AModel, s::NamedTuple)
    vars = variables(model.controls)
    vals = tuple( (model.calibration[k] for k in vars )...)
    return NamedTuple{vars}(vals)
end

# function initial_guess(model::ADModel, m::SVector, s::SVector, p)

#     m = LVectorLike(model.calibration.m, m)  # this does not keep the original type
#     s = LVectorLike(model.calibration.s, s)
    
#     x = initial_guess(model, m, s, p)
    
#     return SVector(x...)
    
# end

function initial_guess(model::AModel, s::SVector)

    ss = NamedTuple{variables(model.states)}(s)
    xx = initial_guess(model, ss)

     return SVector(xx...)

end



function initial_guess(model::AModel, s::QP)
    
    initial_guess(model, s.val)

end


function initial_guess(dmodel::ADModel)
    GVector(
        dmodel.grid,
        [(SVector( initial_guess(dmodel.model, s)... ) for s in Dolo.enum(dmodel.grid))...]
    )
end
