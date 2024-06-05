module DoloCUDAExt

    using CUDA
    using Dolo
# should it be merged with the general definition?
    
    import CUDA: CuArray
    
    Dolo.distance(x::GVector{G, A}, y::GVector{G,A}) where G where A<:CuArray = Base.mapreduce(u->maximum(abs.(u)), max, x.data-y.data)


end
