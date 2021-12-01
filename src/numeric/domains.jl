
### Domains


abstract type AbstractDomain{d} end


struct EmptyDomain{d} <: AbstractDomain{d}
    states::Vector{Symbol}
end

struct CartesianDomain{d} <: AbstractDomain{d}
    states
    min::SVector{d,Float64}
    max::SVector{d,Float64}
end

CartesianDomain(states, min::Vector, max::Vector) = CartesianDomain{length(min)}(
    states,
    SVector(min...),
    SVector(max...)
)

function Base.show(io::IO, cd::CartesianDomain{d}) where d
    rgs = [((cd.min[i]), (cd.max[i])) for i=1:d]
    print(io, replace(join(rgs, "×"), "Inf"=>"∞"))
end

function CartesianDomain(min, max)
    d = length(min)
    states = [string("x", i) for i =1:d]
    return CartesianDomain(states, min, max)
end


struct DiscreteDomain{d} <: AbstractDomain{d}
    points::Vector{Point{d}}
end

struct ProductDomain{d, T,S} <: AbstractDomain{d}
    exo::T
    endo::S
end


function Base.show(io::IO, pd::ProductDomain{d}) where d
    print(io, join( [], "×"))
end

ProductDomain(d1, d2) = ProductDomain{ndims(d1)+ndims(d2), typeof(d1), typeof(d2)}(d1,d2)


function Base.show(io::IO, pd::ProductDomain) where d
    print(io, string(pd.exo, "⊗", pd.endo))
end

⊗(dom1::AbstractDomain{d1}, dom2::AbstractDomain) where d1 where d2 = ProductDomain(dom1, dom2)

struct PDomain{d, T<:Tuple} <: Grid{d}
    domains::T
end

function Base.show(io::IO, pd::PDomain) where d
    print(io, join( [pd.domains], "×"))
end


ndims(t::AbstractDomain{d}) where d = d

ndims(dom::CartesianDomain) = length(dom.min)

discretize(dom::DiscreteDomain{d}) where d = UnstructuredGrid(dom.points)

function discretize(dom::CartesianDomain; n=Union{Int, Vector{Int}})
    if typeof(n)<:Int
        nv = fill(n, ndims(dom))
    else
        nv = n
    end
    min = dom.min
    max = dom.max
    return CartesianGrid(min, max, n)
end
