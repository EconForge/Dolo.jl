
###
# Object to represent a multivalued policy/decision_rule which depends on continuous states only
###

type DecisionRule
    a:: Array{Float64, 1}
    b:: Array{Float64, 1}
    dims:: Array{Int64, 1}
    values:: Array{Float64, 2}
    interpolants:: Any
    function DecisionRule(a, b, dims)
        d = length(dims)
        values = zeros(0,0)
        interpolant = Union{}
        return new(a, b, dims, values, interpolant)
    end
end

function set_values(dr::DecisionRule, values)
    dims = dr.dims
    n_x = size(values)[end]
    d = length(dr.dims)
    knots = [linspace(dr.a[i], dr.b[i], dr.dims[i]) for i=1:d]
    dr.interpolants = [ scale(interpolate((reshape(slice(values,:,i_x),dims...)), BSpline(Linear()), OnGrid()), knots...) for i_x in 1:n_x  ]
end

function evaluate(dr::DecisionRule, s::AbstractArray{Float64,1})
    n_x = length(dr.interpolants)
    return Float64[dr.interpolants[j][s...] for j=1:n_x ]
end

function evaluate(dr::DecisionRule, s::AbstractArray{Float64,2})
    return cat(1, [evaluate(dr, slice(s,n,:))' for n=1:size(s,1)]... )
end

typealias Policy DecisionRule

type ConstantPolicy
    x:: Array{Float64,1}
end

evaluate(cp::ConstantPolicy, s::AbstractArray{Float64,1}) = cp.x
evaluate(cp::ConstantPolicy, s::AbstractArray{Float64,2}) = repmat(cp.x',size(s,1),1)

###
# Object to represent a multivalued decision rule/policy which depends on a discrete state and other continuous states
###

type MixedDecisionRule
    n_mc:: Int64
    a:: Array{Float64, 1}
    b:: Array{Float64, 1}
    dims:: Array{Int64, 1}
    values:: Array{Float64, 3}
    interpolants:: Any
    function MixedDecisionRule(n_mc, a, b, dims)
        d = length(dims)
        values = zeros(0,0,0)
        interpolant = Union{}
        return new(n_mc, a, b, dims, values, interpolant)
    end
end

function set_values(mdr::MixedDecisionRule, values)
    n_mc = mdr.n_mc
    dims = mdr.dims
    n_x = size(values)[end]
    d = length(mdr.dims)
    knots = [linspace(mdr.a[i], mdr.b[i], mdr.dims[i]) for i=1:d]
    mdr.interpolants = [ scale(interpolate((reshape(slice(values,i_mc,:,i_x),dims...)), BSpline(Linear()), OnGrid()), knots...) for i_mc in 1:n_mc, i_x in 1:n_x  ]
end

function evaluate(mdr::MixedDecisionRule, i::Int64, s::AbstractArray{Float64,1})
    n_x = size(mdr.interpolants,2)
    return Float64[mdr.interpolants[i,j][s...] for j=1:n_x ]
end

function evaluate(mdr::MixedDecisionRule, i::Int64, s::AbstractArray{Float64,2})
    return cat(1, [evaluate(mdr, i, slice(s,n,:))' for n=1:size(s,1)]... )
end

typealias MixedPolicy MixedDecisionRule


type MixedConstantPolicy
    x:: Array{Float64,1}
end

evaluate(mcp::MixedConstantPolicy, i, s::AbstractArray{Float64,1}) = mcp.x
evaluate(mcp::MixedConstantPolicy, i, s::AbstractArray{Float64,2}) = repmat(mcp.x',size(s,1),1)
