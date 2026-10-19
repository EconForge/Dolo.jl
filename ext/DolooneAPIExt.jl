module DolooneAPIExt

# Loaded automatically when oneAPI is loaded together with Dolo (see DoloCUDAExt).

using oneAPI: oneArray
import Dolo: GVector, distance
import Dolo.splines: prefilter!

# the generic method iterates with `zip`, which would require scalar indexing
distance(x::GVector{G,<:oneArray}, y::GVector{G,<:oneArray}) where G = mapreduce(distance, max, x.data, y.data)

prefilter!(data::oneArray) = prefilter!(data, Val(:KA))

end
