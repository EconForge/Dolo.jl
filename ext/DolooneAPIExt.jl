module DolooneAPIExt

# Loaded automatically when oneAPI is loaded together with Dolo (see DoloCUDAExt).

using oneAPI: oneArray
import Dolo: GVector, distance

# the generic method iterates with `zip`, which would require scalar indexing
distance(x::GVector{G,<:oneArray}, y::GVector{G,<:oneArray}) where G = mapreduce(distance, max, x.data, y.data)

end
