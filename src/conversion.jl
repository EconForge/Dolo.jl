
function convert_precision(T, model::Dolo.YModel)

    # convert calibration
    calibration = NamedTuple( ((a,T(b)) for (a,b) in pairs(model.calibration) ) )

    # convert exogenous shock
    fun = u->T(u)
    P = fun.(model.exogenous.P)
    Q = SVector((fun.(e) for e in model.exogenous.Q)...)

    vars = Dolo.variables(model.exogenous)
    exogenous = Dolo.MarkovChain(vars, P, Q)

    # convert states and controls spaces
    states = convert_precision(T, model.states)
    controls = convert_precision(T, model.controls)


    Dolo.YModel(
        names,
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
