using StaticArrays: SMatrix

function convert_precision(T, exogenous::Dolo.MarkovChain)

    fun = u->T(u)
    P = fun.(exogenous.P)
    Q = SVector((fun.(e) for e in exogenous.Q)...)
    vars = Dolo.variables(exogenous)
    Dolo.MarkovChain(vars, P, Q)

end

function convert_precision(T, var::Dolo.VAR1)

    d = size(var.Σ,1)

    vars = Dolo.variables(var)

    ρ = convert(T, var.ρ)
    Σ = SMatrix{d,d,T,d*d}(var.Σ)
    
    Dolo.VAR1(vars,ρ,Σ)

end

function convert_precision(T, model::Dolo.YModel)

    # convert calibration
    calibration = NamedTuple( ((a,T(b)) for (a,b) in pairs(model.calibration) ) )

    # convert exogenous shock

    exogenous = convert_precision(T, model.exogenous)

    # convert states and controls spaces
    states = convert_precision(T, model.states)
    controls = convert_precision(T, model.controls)


    Dolo.YModel(
        Dolo.name(model),
        states,
        controls,
        exogenous, 
        calibration
    )
end

function convert_precision(Tnew, cspace::Dolo.CartesianSpace{d, vars, Told}) where Told where d where vars
    Dolo.CartesianSpace{d, vars, Tnew}(cspace.min, cspace.max)
end

function convert_precision(Tnew, gspace::Dolo.GridSpace{N,d,dims,Tf}) where N where d where dims where Tf
    Dolo.GridSpace{N,d,dims,Tnew}(gspace.points)
end

function convert_precision(Tnew, ps::Dolo.ProductSpace)
    Dolo.ProductSpace((convert_precision(Tnew, e) for e in ps.spaces)...)
end


function hypeof(model)
    typ = typeof(model)
    otyp = eltype(model)
    ntyp = Float32
    s = replace(string(typ),string(otyp)=>"Float32")
    gtyp = ( Meta.parse(s) )
    return Union{eval(gtyp), typ}
end


getprecision(model::Dolo.YModel) = typeof(model.calibration[1])
getprecision(dmodel::Dolo.DYModel) = getprecision(dmodel.model)


## Discretized models

convert_precision(T, g::CGrid{d}) where d = CGrid{d,T}(map(r -> (T(r[1]), T(r[2]), r[3]), g.ranges))
convert_precision(T, g::SGrid{n,d}) where n where d = SGrid{n,d,T}(SVector{n}(map(p -> SVector{d,T}(p), g.points)))
convert_precision(T, g::PGrid) = PGrid(convert_precision(T, g.grids[1]), convert_precision(T, g.grids[2]))

# discretized iid shocks: nodes and weights
convert_precision(T, d::NamedTuple{(:x,:w)}) = (; x = map(e -> T.(e), d.x), w = T.(d.w))

"Copy of the discretized model `dmodel` in precision `T` (e.g. Float32)."
convert_precision(T, dmodel::DYModel) = DYModel(
    convert_precision(T, dmodel.model),
    convert_precision(T, dmodel.grid),
    convert_precision(T, dmodel.dproc)
)
