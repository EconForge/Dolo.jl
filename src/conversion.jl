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
