function sample(weights::SVector{d,Float64}) where d
    x = rand()
    findfirst(u-> (u>=x),cumsum(weights))
end

function sample(weights::NTuple{d,Float64}) where d
    sample(SVector(weights...))
end
