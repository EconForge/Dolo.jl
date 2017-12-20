###
# special struct to hold decision rules
###

# import Main
struct DRStore{T}
    data::Tuple{Vararg{Vector{T}}}
    flat::Vector{T}
end

Main.length(ds::DRStore) = length(ds.data)
Main.size(ds::DRStore, x...) = size(ds.data, x...)
Main.getindex(ds::DRStore,x...) = getindex(ds.data,x...)
Main.setindex!(ds::DRStore,x...) = setindex!(ds.data,x...)
Main.abs(ds::DRStore{T}) where T = abs(ds.flat)

Main.maxabs(ds::DRStore{T}) where T = maximum(e->norm(e,Inf), ds.flat)

Main.norm(ds::DRStore{T}) where T = maximum(e->norm(e,Inf), ds.flat)
# Main.norm(ds::DRStore{T},k::Int) where T = norm( (norm(e,k) for e in ds.flat), k)
distance(ds1::DRStore{T}, ds2::DRStore{T}) where T = maximum( k->norm(k[1]-k[2],Inf), zip(ds1.flat,ds2.flat) )
# distance(ds1::DRStore{T}, ds2::DRStore{T}, k::Int) where T = norm( (norm(e[1]-e[2],k) for e in zip(ds1.flat,ds2.flat)), k )

total(ds::DRStore{T}) where T = sum(length(e) for e in ds.data)

function DRStore(::Type{T}, dims::AbstractVector{Int}) where T
    L = sum(dims)
    fD = zeros(T, L)
    ptr = pointer(fD)
    sz = sizeof(T)
    offsets = [[0]; cumsum(dims)[1:end-1]]*sz
    D = [unsafe_wrap(Array{T}, ptr+offsets[i], dims[i]) for i=1:length(dims)]
    return DRStore(tuple(D...), fD)
end
function DRStore(data::AbstractVector{S}) where S<:AbstractVector{T} where T
    # data is always copied in this case
    # as we cannot guarantee data is contiguous
    dims = [length(d) for d in data]
    ds = DRStore(T, dims)
    for i=1:length(data)
        ds.data[i][:] = data[i]
    end
    ds
end

# this one reuses data in place
# should it be called DRStore! ?
function DRStore(fD::Vector{T},dims::AbstractVector{Int}) where T
    L = sum(dims)
    ptr = pointer(fD)
    sz = sizeof(T)
    offsets = [[0]; cumsum(dims)[1:end-1]]*sz
    D = [unsafe_wrap(Array{T}, ptr+offsets[i], dims[i]) for i=1:length(dims)]
    return DRStore(tuple(D...), fD)
end

DRStore! = DRStore

function Main.copy(ds::DRStore{T}) where T
    dims = [length(e) for e in ds.data]
    fD = deepcopy(ds.flat)
    ds1 = DRStore!(fD,dims)
    return ds1
end

Main.copy!(ds0::DRStore{T}, ds1::DRStore{T}) where T = copy!(ds0.flat, ds1.flat)


function Main.clamp!(ds::DRStore{T}, ds_lb::DRStore{T}, ds_ub::DRStore{T}) where T
    N = length(ds.flat)
    for n=1:N
        ds.flat[n] = clamp.(ds.flat[n], ds_lb.flat[n], ds_ub.flat[n])
    end
end
