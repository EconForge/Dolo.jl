struct YModel{N,C,A,B,D,S} <: AModel{C}
    states::A         # must be exo \times endo
    controls::B
    exogenous::C
    calibration::D
    source::S
end 

YModel(N,A,B,C,D) = let
    YModel{N,typeof(C),typeof(A),typeof(B),typeof(D),Nothing}(A,B,C,D,nothing)
end
YModel(N,A,B,C,D,S) = YModel{N,typeof(C),typeof(A),typeof(B),typeof(D),typeof(S)}(A,B,C,D,S)

name(::YModel{N,C,A,B,D}) where C where A where B where D where N = N

bounds(model::YModel, s) = model.controls


# TODO: check somewhere that the type of
#  states/controls/exogenous
# is the same as calibration

eltype(model::YModel) = eltype(model.calibration)

function recalibrate(model::YModel; kwargs...)
    calib = merge(model.calibration, kwargs)
    YModel(model.states, model.controls, model.exogenous, calib, model.source)
end

get_states(model::YModel) = variables(model.states)
get_controls(model::YModel) = variables(model.controls)
# get_endo_states(model::Dolo)
# get_exo_states(model::Dolo)

# discretize(cc::CartesianSpace; n=10) = CGrid( tuple(( (cc.min[i],cc.max[i], n) for i=1:length(cc.min))...) )


import Term: Panel, tprint, tprintln

function Base.show(io::IO, m::YModel) 
    # println(Panel("this is {red}RED{/red}"; fit=true))
    
    exovars = [e for e in Dolo.variables(m.exogenous)]
    exotype = typeof(m.exogenous).name.name
    hcontrols = join(get_controls(m), ", ")
    hstates = join([
        (e in exovars ? "{red}$e{/red}" : e)
        for e in get_states(m)
    ], ", ")

    txt = """
    Model: {blue}$(name(m)){/blue}
    * states: $hstates
    * controls: $hcontrols
    * exogenous: $exotype({red}$(join(exovars,",")){/red})"""

    tprintln(Panel(txt; fit=true))
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

# TODO: check whether true
eltype(dm::DYModel) = eltype(dm.model)


name(dm::DYModel) = name(dm.model)

bounds(dmodel::DYModel, s) = bounds(dmodel.model, s)

function discretize(model::AModel{<:MvNormal}, d=Dict())
    exo = get(d, :exo, Dict())
    endo = get(d, :endo, Dict())

    dist = discretize(model.exogenous, exo)
    grid = discretize(model.states, endo)
    return DYModel(model, grid, dist)
end

function discretize(model::AModel{<:VAR1}, d=Dict())

    exo = get(d, :exo, Dict())
    endo = get(d, :endo, Dict())
    dvar = discretize(model.exogenous, exo)
    d = size(model.exogenous.Σ,1)

    # number of endogenous states
    n_s = length(Dolo.variables(model.states)) - d
    
    exo_grid = SGrid(dvar.Q)
    # TODO: simplify that call
    Tf = eltype(model)
    endo_space = CartesianSpace{n_s, Dolo.variables(model.states)[d+1:end],Tf}(
        model.states.endo.min,
        model.states.endo.max
    )
    # return endo_space
    endo_grid = discretize(endo_space, endo)
    grid = exo_grid × endo_grid
    return DYModel(model, grid, dvar)
end

function discretize(model::AModel{<:MarkovChain}, d=Dict())
    # exo = get(d, :exo, Dict())
    endo = get(d, :endo, Dict())
    dvar = discretize(model.exogenous)
    exo_grid = SGrid(dvar.Q)
    endo_grid = discretize(model.states.spaces[2], endo)
    grid = exo_grid × endo_grid
    return Dolo.DYModel(model, grid, dvar)
end
