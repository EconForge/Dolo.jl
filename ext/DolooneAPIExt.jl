module DolooneAPIExt
    
    println("Loading oneAPI")
    using oneAPI
    import oneAPI: oneArray

    using Dolo: distance

    Dolo.distance(x::GVector{G, A}, y::GVector{G,A}) where G where A<:oneArray = Base.mapreduce(u->maximum(abs.(u)), max, x.data-y.data)

    import Dolo.splines: prefilter!

    prefilter!(data::oneArray) = prefilter!(data, Val(:KA))


end
