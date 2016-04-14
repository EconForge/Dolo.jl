using QuantEcon: gridmake

#####
# helper functions/types (to be moved/removed)
#####

mlinspace(a,b,dims) = gridmake([linspace(a[i],b[i],dims[i]) for i=1:length(dims)]...)


type MarkovChain
    P::Array{Float64, 2}
    Q::Array{Float64, 2}
end

type ApproximationSpace
  a:: Array{Float64, 1}
  b:: Array{Float64, 1}
  dims:: Array{Int64, 1}
end
