abstract type Space{d} end

struct CartesianSpace{d,dims}
    # names::NTuple{d, Symbol}
    min::NTuple{d, Float64}
    max::NTuple{d, Float64}
end

const CSpace = CartesianSpace

CartesianSpace(a::Tuple{Float64}, b::Tuple{Float64}) = CartesianSpace{length(a), (:x,)}(a,b)
CartesianSpace(a::Tuple{Float64, Float64}, b::Tuple{Float64, Float64}) = CartesianSpace{length(a), (:x_1, :x_2)}(a,b)

function CartesianSpace(;kwargs...)
    names = tuple(keys(kwargs)...)
    a = tuple((v[1] for v in values(kwargs))...)
    b = tuple((v[2] for v in values(kwargs))...)
    d = length(names)
    return CartesianSpace{d, names}(a,b)
end

getindex(cs::CartesianSpace{d}, ind::SVector{d, Float64}) where d = ind

import Base: in

function draw(cs::CartesianSpace) 
    loc = SVector(( (cs.min[i] + rand()*(cs.max[i]-cs.min[i])) for i=1:length(cs.min) )...)
    val = loc
    QP(loc,val)
    # (;loc, val)
end

Base.in(e::SVector, cs) = all( ( (e[i]<=cs.max[i])&(e[i]>=cs.min[i]) for i=1:length(cs.min) ) )

# TODO: why is this not working?
# dims(dom::CartesianSpace{d,dims}) where d where dims = dims


ndims(dom::CartesianSpace{d, dims}) where d where dims = d
variables(dom::CartesianSpace{d,t}) where d where t = t
dims(dom::CartesianSpace) = variables(dom)

struct GridSpace{N,d,dims}
    points::SVector{N,SVector{d,Float64}}
end

const GSpace = GridSpace

GridSpace(v::SVector{N, SVector{d, Float64}}) where d where N = GridSpace{length(v), d, (:i_,)}(SVector(v...))
GridSpace(v::Vector{SVector{d, Float64}}) where d = GridSpace{length(v), d, (:i_,)}(SVector(v...))
GridSpace(names, v::SVector{k, SVector{d, Float64}}) where k where d = GridSpace{length(v), d, names}(v)

getindex(gs::GridSpace, i::Int64) = gs.points[i]

function draw(g::GridSpace) 
    i =  rand(1:length(g.points))   # loc
    v = g.points[i]                 # val
    QP(i,v)
    # (;loc=i,val=v)
end

ndims(gd::GridSpace{N,d,dims}) where N where d where dims = d
ddims(gd::GridSpace{N,d,dims}) where N where d where dims = dims
dims(gd::GridSpace) = ddims(gd)


struct ProductSpace{A,B}
    spaces::Tuple{A,B}
end

ProductSpace(A,B) = ProductSpace((A,B))



function draw(p::ProductSpace)
    a = rand(p.spaces[1])
    b = rand(p.spaces[2])
    QP(
        (a.loc,b.loc),
        SVector(a.val...,b.val...)
    )
    # (;loc=(a.loc,b.loc),val=SVector(a.val...,b.val...))
end

getindex(ps::ProductSpace, ind) = SVector(getindex(ps.spaces[1], ind[1])..., getindex(ps.spaces[2], ind[2])...)

variables(p::ProductSpace) = dims(p)

import LinearAlgebra: cross
cross(A::DA, B::DB) where DA<:GridSpace where DB<:CartesianSpace = ProductSpace{DA, DB}((A,B))

# (×)(A::DA, B::DB) where DA<:GridSpace where DB<:CartesianSpace = cross(A,B)

ndims(p::P) where P<:ProductSpace = ndims(p.spaces[1]) + ndims(p.spaces[2])
dims(p::P) where P<:ProductSpace = tuple(dims(p.spaces[1])..., dims(p.spaces[2])...)




# construct QP poitns

function dropnames(namedtuple::NamedTuple, names::Tuple{Vararg{Symbol}}) 
    keepnames = Base.diff_names(Base._nt_names(namedtuple), names)
   return NamedTuple{keepnames}(namedtuple)
end


QP(space::CartesianSpace{d}; values...) where d = let
    s_ =zero(SVector{d,Float64})*NaN
    s0 = QP(s_,s_)
    QP(space, s0; values...)
end

function QP(space::CartesianSpace, s0::QP; values...)

    vars = Dolo.variables(space)

    @assert (keys(values) ⊆ vars)

    nloc = NamedTuple{vars}(s0.loc)
    mloc = merge(nloc, values)
    loc = SVector(mloc...)

    return QP(loc, loc)

end


QP(space::GridSpace{d}; i_=1) where d = QP(i_,space[i_])



function QP(space::ProductSpace{<:GridSpace{d1},<:CartesianSpace{d2}}; values...) where d1 where d2
    i0 = get(values, :i_, 1)
    m_ = space.spaces[1][i0]
    s_ =zero(SVector{d2,Float64})*NaN
    s0 = QP((i0,s_),SVector(m_...,s_...))
    QP(space, s0; values...)
end

function QP(space::ProductSpace{<:GridSpace{d1},<:CartesianSpace{d2}}, s0::QP; values...) where d1 where d2


    vars = Dolo.variables(space.spaces[2])

    @assert (keys(values) ⊆ tuple(:i_,vars...))
    
    
    if :i_ in keys(values)
        i_ = values[:i_]
        vals = NamedTuple(k=>v for (k,v) in values if k!=:i_)
    else
        i_ = s0.loc[1]
        vals = values
    end

    nloc = NamedTuple{vars}(s0.loc[2])

    mloc = merge(nloc, vals)

    loc = (i_, SVector(mloc...))
    vv = space[loc]

    return QP(loc, vv)

end


### discretization



const DEFAULT_GRID_NPOINTS = 20

# discretize(space::CSpace; kwargs...)  = discretize(space, DEFAULT_GRID_NPOINTS)
# discretize(space::CSpace{d}, n::Int) where d = discretize(space, tuple( (n for i=1:d)...))

function discretize(space::CSpace{d}, n::NTuple{d, <:Int}) where d
    CGrid(
        # variables(space),
        tuple( ((space.min[i], space.max[i], n[i]) for i=1:d)... )
    )
end

discretize(space::CSpace, n::Vector) = discretize(space, tuple(n...))

function discretize(space::CSpace, n=DEFAULT_GRID_NPOINTS; nvals...)
    nn = tuple( ( get(nvals, v, n)   for v in variables(space) )...)
    discretize(space, nn)
end

# compat call
function discretize(space::CSpace, d::Dict)
    k = get(d, :n, DEFAULT_GRID_NPOINTS)
    discretize(space, k)
end

#### 
#### Discretize Grid Spaces
####

function discretize(space::GSpace, args...; kwargs...)
    SGrid(space.points)
end

#### 
#### Discretize Product Spaces
####

###
### CSpace × CSpace

function discretize(ps::ProductSpace{<:CSpace{d1}, <:CSpace{d2}}, n::NTuple{d1d2, <:Int}) where d1 where d2 where d1d2
    @assert d1d2==d1+d2
    s1,s2 = ps.spaces
    n1 = tuple((n[i] for i=1:d1)...)
    n2 = tuple((n[i] for i=(d1+1):(d1+d2))...)
    g1 = discretize(s1,n1)
    g2 = discretize(s2,n2)
    ProductGrid(g1, g2)
end

function discretize(ps::ProductSpace{<:CSpace{d1}, <:CSpace{d2}}, n::Int) where d1 where d2
    nn = tuple( (n for i=1:(d1+d2))...)
    discretize(ps, nn)
end

function discretize(ps::ProductSpace; kwargs...)
    ProductGrid(discretize(ps.spaces[1]; kwargs...), discretize(ps.spaces[2]; kwargs...))
end

