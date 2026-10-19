# abstract model
abstract type AbstractModel{Exo} end
abstract type AbstractDModel end

const AModel = AbstractModel
const ADModel = AbstractDModel

macro akin(model)
    typ = typeof(eval(model))
    ntyp = Meta.parse(repr(typ))
    ntyp.args = ntyp.args[1:end-1]
    return ntyp
end


include("dolo_model.jl")

include("model_extensions.jl")

# TODO: restore initial_guess mechanism

# default implementation
function initial_guess(model::AModel, s::NamedTuple)
    vars = variables(model.controls)
    vals = tuple( (model.calibration[k] for k in vars )...)
    return NamedTuple{vars}(vals)
end


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
