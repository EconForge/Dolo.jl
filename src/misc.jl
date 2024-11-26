function sample(weights::SVector{d,Tf}) where d where Tf
    x = Tf(rand()) # here rand returns a Float64, check whether there is a version for any type
    findfirst(u-> (u>=x),cumsum(weights))
end

function sample(weights::NTuple{d,Tf}) where d where Tf
    sample(SVector(weights...))
end
