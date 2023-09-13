function tabulate(model::YModel, dr, key::Symbol, K=100; kwargs...)
    s0 = calibrated(QP, model, :states)
    states = model.states

    s0 = QP(states, s0; kwargs...)
    if typeof(states) <: CartesianSpace
        vars = variables(states)
        i = findfirst(u->u==key, vars)
        a = states.min[i] 
        b = states.max[i] 
    elseif typeof(states) <: ProductSpace
        vars = variables(states.spaces[2])
        i = findfirst(u->u==key, vars)
        a = states.spaces[2].min[i] 
        b = states.spaces[2].max[i] 
    end
    vals = []
    rr = range(a,b;length=K)
    for v in rr
        dd = Dict(key=>v)
        s = QP(states, s0; dd...)
        x = dr(s)
        push!(vals, (s,x))
    end
    H = [SVector(e[1].val..., e[2]...) for e in vals]
    M = vcat( (e' for e in H)... )
    # vars = cat(variables(model.states)..., variables(model.controls)...; dims=1)

    vv = variables(dr)
    if typeof(vv)<:Symbol
        vvv = (vv,)
    else
        vvv = vv
    end

    vars = cat(variables(model.states)..., vvv...; dims=1)
    dd = Dict(
        key => rr,
        :V => vars
    )

    AxisArray(M; dd...)
end

function tabulate(dr::Union{Dolo.Fun, Dolo.DFun}, key::Symbol, s0::QP, K=100; kwargs...)
    states = dr.domain
    # s0 = QP(states; kwargs...)
    if typeof(states) <: CartesianSpace
        vars = variables(states)
        i = findfirst(u->u==key, vars)
        a = states.min[i] 
        b = states.max[i] 
    elseif typeof(states) <: ProductSpace
        vars = variables(states.spaces[2])
        i = findfirst(u->u==key, vars)
        a = states.min[i] 
        b = states.max[i] 
    end
    vals = []
    rr = range(a,b;length=K)
    for v in rr
        dd = Dict(key=>v)
        s = QP(states, s0; dd...)
        x = dr(s)
        push!(vals, (s,x))
    end
    H = [SVector(e[1].loc..., e[2]...) for e in vals]
    M = vcat( (e' for e in H)... )
    vars = range(1, size(M,2))
    dd = Dict(
        key => rr,
        :V => vars
    )
    AxisArray(M; dd...)
    # TODO:  to label columns, the DFun, needs to know output space
end



function calibrated(::Type{QP}, model::YModel, group)
    @assert group==:states
    sv = calibrated(model, :states)
    if typeof(model.states) <: CartesianSpace
        QP(sv, sv)
    elseif typeof(model.states) <: ProductSpace
        d = length(variables(model.exogenous))
        vals = SVector( (sv[i] for i=(d+1:length(sv)))...)
        loc = (1,vals)
        QP(loc, sv)
        # or QP(loc, model.states[loc]) ???
    end
end

function calibrated(model::YModel, group)
    if group==:states
        vars = variables(model.states)
    elseif group==:controls
        vars = variables(model.controls)
    elseif group==:exogenous
        vars = variables(model.exogenous)
    elseif group==:parameters
       return SVector(
            (v for (k,v) in model.calibration if 
                !(k in union(variables(model.states), variables(model.controls), variables(model.exogenous)))
            )...
       )
    else
        error("Unknown group '$(group)'.")
    end
    SVector( (get(model.calibration,v,NaN) for v in vars)... )
end

function calibrated(model::YModel)
    states = calibrated(model, :states)
    controls = calibrated(model, :controls)
    exogenous = calibrated(model, :exogenous)
    (;states, controls, exogenous)
end

QP(space::Space, loc::QP) = loc
QP(space::Space, loc::SVector) = (loc, space[loc])
QP(space::Space, loc::Vector) = QP(space, SVector(loc...))
get_QP(space, loc) = QP(space, loc)

using AxisArrays

function simulate(model::YModel, φ; T=100, kwargs...)

    s0 = calibrated(QP, model, :states)
    states = model.states
    s0 = QP(states, s0; kwargs...)
    simulate(model, φ, s0; T=T)

end


function simulate(model::YModel, φ, s0; T=100)
    # s0 = calibrated_values(model)
    s = Dolo.get_QP(model.states, s0)
    x = φ(s)
    sim = [(s,x)]
    for t=1:T
        sn = transition(model, s, x)
        xn = φ(sn)
        push!(sim, (sn, xn))
        s = sn
        x = xn
    end
    arr = [SVector(e[1].val..., e[2]...) for e in sim]
    vars = tuple(variables(model.states)..., variables(model.controls)...)
    return AxisArray(
        hcat((e for e in arr)...)',
        T = 0:T,
        V = [vars...],
    )
end

function simulate(dmodel::DYModel{M,G,P}, φ, s0; T=100) where M where G where P<:MarkovChain

    # s0 = calibrated_values(model)
    # s = Dolo.get_QP(dmodel.model.states, s0)
    s = s0
    # return "HI"
    x = φ(s)
    sim = [(s,x)]
    for t=1:T

        # super inefficient
        options = tuple( τ(dmodel, s, x)... )
        weights = tuple( (e[1] for e in options) ...)
        i = sample(weights)
        sn = options[i][2]
        xn = φ(sn)
        push!(sim, (sn, xn))
        s = sn
        x = xn
    end
    arr = [SVector(e[1].val..., e[2]...) for e in sim]
    vars = tuple(variables(dmodel.model.states)..., variables(dmodel.model.controls)...)
    return AxisArray(
        hcat((e for e in arr)...)',
        T = 0:T,
        V = [vars...],
    )
end