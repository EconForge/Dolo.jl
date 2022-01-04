function compute_matrix(L::LinearMap)
    M = L.M
    N = L.N
    cols = [L*[k==j ? 1.0 : 0.0 for k=1:N] for j=1:N]
    return cat(cols...; dims=2)
end