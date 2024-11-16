module DoloCUDAExt

    using CUDA
    # using Dolo
# should it be merged with the general definition?
    import Dolo: GVector, distance
    # import CUDA: CuArray
    
    distance(x::GVector{G, A}, y::GVector{G,A}) where G where A<:CuArray = Base.mapreduce(u->maximum(abs.(u)), max, x.data-y.data)
    norm(x::GVector{G, A}) where G where A<:CuArray = Base.mapreduce(u->maximum(abs.(u)), max, x.data)


end
