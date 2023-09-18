@inline function solve_coefficients!(
            bands::AbstractMatrix{Float64},
            bb::AbstractVector{T},
            coefs::AbstractVector{T}
    ) where T

    M = size(coefs,1) - 2

    # Solve interpolating equations
    # First and last rows are different

    bands[1, 2] /= bands[1, 1]
    bands[1, 3] /= bands[1, 1]
    bb[1] /= bands[1, 1]
    bands[1, 1] = 1.0
    bands[2, 2] -= bands[2, 1] * bands[1, 2]
    bands[2, 3] -= bands[2, 1] * bands[1, 3]
    bb[2] -= bands[2, 1] * bb[1]
    bands[1, 1] = 0.0
    bands[2, 3] /= bands[2, 2]
    bb[2] /= bands[2, 2]
    bands[2, 2] = 1.0

    # Now do rows 2 through M+1
    for row = 3:(M+1)
        bands[row, 2] -= bands[row, 1] * bands[row - 1, 3]
        bb[row] -= bands[row, 1] * bb[row - 1]
        bands[row, 3] /= bands[row, 2]
        bb[row] /= bands[row, 2]
        bands[row, 1] = 0.0
        bands[row, 2] = 1.0
    end
    # Do last row
    bands[M + 2, 2] -= bands[M + 2, 1] * bands[M, 3]
    bb[M + 2] -= bands[M + 2, 1] * bb[M]
    bands[M + 2, 3] -= bands[M + 2, 2] * bands[M + 1, 3]
    bb[M + 2] -= bands[M + 2, 2] * bb[M + 1]
    bb[M + 2] /= bands[M + 2, 3]
    bands[M + 2, 3] = 1.0

    coefs[M + 2] = bb[M + 2]
    # Now back substitute up
    # for row in (M+1):2
    for ii in 1:M
        row = M+2 - ii
        coefs[row] = bb[row] - bands[row, 3] * coefs[row + 1]
    end
    # Finish with first row
    coefs[1] = bb[1] - bands[1, 2] * coefs[2] - bands[1, 3] * coefs[3]

end

@inline function fill_bands!(M::Int, bands::AbstractMatrix, bb::AbstractVector{T}, data::AbstractVector{T}) where T

    # boundary conditions
    @inbounds bands[1, 1] = 1.0
    @inbounds bands[1, 2] = -2.0
    @inbounds bands[1, 3] = 1.0
    @inbounds bands[M+2, 1] = 1.0
    @inbounds bands[M+2, 2] = -2.0
    @inbounds bands[M+2, 3] = 1.0

    for k=2:M+1
        @inbounds bands[k,1] = 1.0/6.0
        @inbounds bands[k,2] = 2.0/3.0
        @inbounds bands[k,3] = 1.0/6.0
    end

    @inbounds bb .= data
    @inbounds bb[1] = zero(T)
    @inbounds bb[M+2] = zero(T)

end


# TODO: there are certainly some cases where we don't want MVectors there

function prefilter!(data::AbstractVector{T}) where T

    N = length(data)

    bands = zero(MMatrix{N, 3, Float64, N*3})
    bb = zero(MVector{N, T})

    fill_bands!(N-2, bands, bb, data)
    solve_coefficients!(bands, bb, data)

end


# function prefilter!(data::AbstractVector{T}) where T

#     # TODO: find a way to avoid allocation

#     # data is 1d array of length M
#     # data is 1d array of length M+2
#     M = length(data) - 2
#     bands = zeros(M+2, 3)
    
#     bb = zeros(T, M+2)

#     # fill_bands!(M, bands, bb, data)
#     prefilter!(data, bands, bb)

# end

# function prefilter!(data, bands, bb)

#     M = length(data) - 2
#     fill_bands!(M, bands, bb, data)
#     solve_coefficients!(bands, bb, data)

# end



# function prefilter!(data::Array{T,1}) where T
#     M = size(data,1)-2
#     bands = zeros(M+2,3)
#     bb = zeros(T,M+2)
#     fill_bands!(M, bands, bb, data)
#     prefilter!(data, bands, bb)
# end

function prefilter!(data::AbstractArray{T,2}) where T

    I,J = size(data)

    M = J-2
    for i=1:I
        dat = view(data, i, :)
        prefilter!(dat)
    end

    M = I-2
    for j=1:J
        dat = view(data, :, j)
        prefilter!(dat)
    end
end

function prefilter!(data::Array{T,3}) where T

    I,J,K = size(data)

    M = K-2
    # for i=1:I
    for i=1:I
        for j=1:J
            dat = view(data, i,j, :)
            prefilter!(dat)
        end
    end
    M = J-2
    # Threads.@threads for i=1:I
    for i=1:I
        for k=1:K
            dat = view(data, i,:, k)
            prefilter!(dat)
        end
    end
    M = I-2
    # Threads.@threads for j=1:J
    for j=1:J
        for k=1:K
            dat = view(data, : ,j, k)
            prefilter!(dat)
        end
    end
end


function prefilter!(data::Array{T,4}) where T
    I,J,K,L = size(data)

    M = L-2

    for i=1:I
        for j=1:J
            for k=1:K
                dat = view(data,i,j,k,:)
                prefilter!(dat)
            end
        end
    end
    M = K-2
    for i=1:I
        for j=1:J
            for l=1:L
                dat = view(data,i,j,:,l)
                prefilter!(dat)
            end
        end
    end
    M = J-2
    for i=1:I
        for k=1:K
            for l=1:L
                dat = view(data, i,:, k,l)
                prefilter!(dat)
            end
        end
    end
    M = I-2
    for j=1:J
        for k=1:K
            for l=1:L
                dat = view(data, :, j,k,l)
                prefilter!(dat)
            end
        end
    end
end
