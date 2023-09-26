struct YModel{C,A,B,D,N,S} <: AModel
    states::A         # must be exo \times endo
    controls::B
    exogenous::C
    calibration::D
    source::S
end 

YModel(N,A,B,C,D) = let
    println("Who is calling?")
    YModel{typeof(C),typeof(A),typeof(B),typeof(D),N,Nothing}(A,B,C,D,nothing)
end
YModel(N,A,B,C,D,S) = YModel{typeof(C),typeof(A),typeof(B),typeof(D),N,typeof(S)}(A,B,C,D,S)

name(::YModel{C,A,B,D,N}) where C where A where B where D where N = N

bounds(model::YModel, s) = model.controls

function recalibrate(model::YModel; kwargs...)
    calib = merge(model.calibration, kwargs)
    YModel(model.states, model.controls, model.exogenous, calib, model.source)
end

get_states(model::YModel) = variables(model.states)
get_controls(model::YModel) = variables(model.controls)
# get_endo_states(model::Dolo)
# get_exo_states(model::Dolo)

# discretize(cc::CartesianSpace; n=10) = CGrid( tuple(( (cc.min[i],cc.max[i], n) for i=1:length(cc.min))...) )


function Base.show(io::IO, m::YModel) 
    println("Model")
    println("* name: ", name(m))
    println("* states: ", join(get_states(m), ", "))
    println("* controls: ", join(get_controls(m), ", "))
    println("* exogenous: ", m.exogenous)
end

function Base.show(io::IO, m::ADModel) 
    println("Discretized Model")
    println("* name = ", name(m))
end


struct DYModel{M, G, D} <: ADModel
    model::M
    grid::G
    dproc::D
end

name(dm::DYModel) = name(dm.model)

bounds(dmodel::DYModel, s) = bounds(dmodel.model, s)

function discretize(model::YModel{<:MvNormal}, d=Dict())
    exo = get(d, :exo, Dict())
    endo = get(d, :endo, Dict())

    dist = discretize(model.exogenous, exo)
    grid = discretize(model.states, endo)
    return DYModel(model, grid, dist)
end

function discretize(model::YModel{<:VAR1}, d=Dict())
    exo = get(d, :exo, Dict())
    endo = get(d, :endo, Dict())
    dvar = discretize(model.exogenous, exo)
    d = size(model.exogenous.Σ,1)
    n_s = length(Dolo.variables(model.states)) - d
    
    exo_grid = SGrid(dvar.Q)
    endo_space = CartesianSpace{n_s, Dolo.variables(model.states)[d+1:end]}(
        model.states.min[d+1:end],
        model.states.max[d+1:end]
    )
    # return endo_space
    endo_grid = discretize(endo_space, endo)
    grid = exo_grid × endo_grid
    return DYModel(model, grid, dvar)
end

function discretize(model::YModel{<:MarkovChain}, d=Dict())
    # exo = get(d, :exo, Dict())
    endo = get(d, :endo, Dict())
    dvar = discretize(model.exogenous)
    exo_grid = SGrid(dvar.Q)
    endo_grid = discretize(model.states.spaces[2], endo)
    grid = exo_grid × endo_grid
    return Dolo.DYModel(model, grid, dvar)
end
