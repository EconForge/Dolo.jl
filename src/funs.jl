struct Policy{From, To, Fun}
    from::From
    to::To
    fun::Fun
end

(pol::Policy)(s::QP) = pol.fun(s.val)
(pol::Policy)(s) = pol.fun(s)


using .splines
using .splines: Linear, Cubic, MLinear, MCubic, CubicInterpolator, SplineInterpolator



struct DFun{Dom, Gar, Itp, vars}
    domain::Dom
    values::Gar
    itp::Itp
end

vars(::DFun{D,G,I,vs}) where D where G where I where vs = vs
variables(::DFun{D,G,I,vs}) where D where G where I where vs = vs


DFun(domain, values, itp, vars) = DFun{typeof(domain), typeof(values), typeof(itp), vars}(domain, values, itp)
function DFun(domain, values, itp)
    if eltype(values) <: Number
        vars = :y
    elseif eltype(values) <: SVector # && length(eltype(values)) == 1
        vars = (:y,)
    else
        println(values)
    end
    return DFun(domain, values, itp, vars)
end
# DFun(domain, values, itp) = begin @assert ndims(domain)==1 ; DFun(domain, values, itp, (:y,)) end
# works for cartesian only
function DFun(states, values::GVector{G,V}, vars=nothing; interp_mode=:linear) where V where G<:CGrid

    vv = reshape(values.data, (e[3] for e in values.grid.ranges)...)
    if interp_mode == :linear
        itp = SplineInterpolator(values.grid.ranges; values=vv, k=1)
    elseif interp_mode == :cubic
        itp = SplineInterpolator(values.grid.ranges;  values=vv,k=3)       
    else
        throw("Unkown interpolation mode $(interp_mode)")
    end

    if typeof(vars) <: Nothing
        if eltype(values) <: Number
            vars = :y
        elseif eltype(values) <: SVector # && length(eltype(values)) == 1
            vars = (:y,)
        else
            println(values)
        end
    end

    return DFun(states, values, itp, vars)
end

function DFun(states, values::GVector{G,V}, vars=nothing; interp_mode=:linear) where V where G<:PGrid{G1,G2} where G1<:SGrid where G2<:CGrid
    
    if interp_mode == :linear
        k=1
    elseif interp_mode == :cubic
        k=3
    else
        throw("Unkown interpolation mode $(interp_mode)")
    end

    # TODO: check values.data[i,:]
    sz = (e[3] for e in values.grid.g2.ranges)
    itps = tuple( (SplineInterpolator(values.grid.g2.ranges;  values=reshape(values[i,:], sz...),k=k)  for i=1:length(values.grid.g1)  )...)

    if typeof(vars) <: Nothing
        if eltype(values) <: Number
            vars = :y
        elseif eltype(values) <: SVector # && length(eltype(values)) == 1
            vars = (:y,)
        else
            println(values)
        end
    end

    return DFun(states, values, itps, vars)

end

function fit!(φ::DFun, x::GVector{G}) where G<:PGrid{<:SGrid, <:CGrid}

    # This is only for SGrid x CGrid

    n_m = length(x.grid.g1)
    sz = tuple(n_m, (e[3] for e in x.grid.g2.ranges)...)    
    xx = reshape( view(x.data, :), sz...)
    rr = tuple( (Colon() for i=1:(Base.ndims(xx)-1))... )
    for i=1:length( φ.itp)
        vv = view( xx, i, rr...)
        splines.fit!(φ.itp[i], vv)
    end

end

function fit!(φ::DFun, x::GVector{G}) where G<:CGrid

    splines.fit!(φ.itp, x.data)

end

## PGrid
function (f::DFun{A,B,I,vars})(i::Int, x::SVector{d2, U})  where A where B<:GArray{G,V} where V where I where G<:PGrid{G1,G2} where G1<:SGrid where G2<:CGrid where vars where d2 where U
    f.itp[i](x)
end

function (f::DFun{A,B,I,vars})(x::QP)  where A where B<:GArray{G,V} where V where I where G<:PGrid{G1,G2} where G1<:SGrid where G2<:CGrid where vars
    f(x.loc...)
end

function (f::DFun{A,B,I,vars})(x::Tuple)  where A where B<:GArray{G,V} where V where I where G<:PGrid{G1,G2} where G1<:SGrid where G2<:CGrid where vars
    f(x...)
end

function (f::DFun{A,B,I,vars})(loc::Tuple{Tuple{Int64}, SVector{d2, U}})  where A where B<:GArray{G,V} where V where I where G<:PGrid{G1,G2} where G1<:SGrid where G2<:CGrid where vars where d2 where U
    # TODO: not beautiful
    x_ = loc[2]
    dd1 = ndims(f.values.grid.g1)
    dd2 = ndims(f.values.grid.g2)
    x = SVector((x_[i] for i=(dd1+1):(dd1+dd2))...)
    f.itp[loc[1][1]](x)
end

# CGrid

function (f::DFun{A,B,I,vars})(x::QP)  where A where B<:GArray{G,V} where V where I where G<:CGrid where vars
    return f.itp(x.loc)
end

function (f::DFun{A,B,I,vars})(x::SVector)  where A where B<:GArray{G,V} where V where I where G<:CGrid where vars
    return f.itp(x)
end

# function (f::DFun{A,B,I,vars})(x::Tuple)  where A where B<:GArray{G,V} where V where I where G<:CGrid where vars
#     return f.itp(x)
# end


# Compatibility calls

(f::DFun)(x::Float64) = f(SVector(x))
(f::DFun)(x::Float64, y::Float64) = f(SVector(x,y))
(f::DFun)(x::Vector{SVector{d,Float64}}) where d = [f(e) for e in x]


ndims(df::DFun) = ndims(df.domain)




##
## functions defined by hand
##

struct Fun{Dom, FF}
    domain::Dom
    f::FF
end

Fun(f) = Fun(
    begin n=methods(f)[1].nargs-1; CartesianSpace( tuple( (-Inf for i=1:n)...), tuple( (Inf for i=1:n)...) ) end,
    f
)

function (f::Fun{Dom})(args...) where Dom
    return f.f(args...)
end

# function ddFun(dom, gar)
#     d = ndims(dom)
#     dt = typeof(dom)
#     gt = typeof(gar)
#     DFun{d, dt, gt}(dom, gar)
# end

# function discretize(space::D2; n=ntuple(u->15,d))  where  D2<:CartesianSpace{d} where d
    
#     if false in (isfinite.(space.max) .& isfinite.(space.min))
#         # TODO : improve
#         throw(DomainError("Impossible to discretize space with infinite bounds."))
#     end

#     grid = CGrid(
#         tuple( 
#             ( (space.min[i], space.max[i], n[i]) for i=1:ndims(space) )...
#         )
#     )

#     return grid

# end


# function discretize(dom::ProductSpace{D1,D2}; n=ntuple(u->15,ndims(dom))) where D1<:GridSpace where  D2<:CartesianSpace
#     grid_a = SGrid(dom.spaces[1].points)
#     g_b = dom.spaces[2]
#     grid_b = CGrid(
#         tuple( 
#             ( (g_b.min[i], g_b.max[i], n[i]) for i=1:ndims(g_b) )...
#         )
#     )
#     return grid_a × grid_b
# end
