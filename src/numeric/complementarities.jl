function φmin(u::Point{d},v::Point{d}) where d
    sq = sqrt.(u.^2+v.^2)
    # p = (u+v-sq)/2
    # J_u = (1 .- u./sq)/2
    # J_v = (1 .- v./sq)/2
    p =    SVector{d}( (v[i]<Inf ? (u[i]+v[i]-sq[i])/2 : u[i]) for i=1:d )
    J_u =  SVector{d}( (v[i]<Inf ? (1.0 - u[i]./sq[i])/2 : 1) for i=1:d )
    J_v =  SVector{d}( (v[i]<Inf ? (1.0 - v[i]./sq[i])/2 : 0) for i=1:d )

    return p, SDiagonal(J_u), SDiagonal(J_v)
end

function φmax(u::Point{d},v::Point{d}) where d
    sq = sqrt.(u.^2+v.^2)
    # p = (u+v+sq)/2
    # J_u = (1 .+ u./sq)/2
    # J_v = (1 .+ v./sq)/2

    p =    SVector{d}( (v[i]>-Inf ? (u[i]+v[i]+sq[i])/2 : u[i]) for i=1:d )
    J_u =  SVector{d}( (v[i]>-Inf ? (1.0 + u[i]./sq[i])/2 : 1) for i=1:d )
    J_v =  SVector{d}( (v[i]>-Inf ? (1.0 + v[i]./sq[i])/2 : 0) for i=1:d )

    return p, SDiagonal(J_u), SDiagonal(J_v)
end

function PhiPhi(f::Point{d}, x::Point{d}, a::Point{d}, b::Point{d}) where d
    y, y_f, y_x = φmin(f,x-a)

    z, z_y, z_x = φmax(y,x-b)
    return z, (z_y*y_f), (z_y*y_x + z_x)
end

function PhiPhi0(f::Point{d}, x::Point{d}, a::Point{d}, b::Point{d}) where d
    y = φmin(f,x-a)[1]
    z = φmax(y,x-b)[1]
    return z
end


function PhiPhi(f::Point{d}, D::SMatrix{d,d,Float64,q}, x::Point{d}, a::Point{d}, b::Point{d}) where d where q
    z, z_f, z_x = PhiPhi(f,x,a,b)
    z, z_f*D + z_x
end



function PhiPhi!(F::Vector{Vector{Point{d}}},
                 X::Vector{Vector{Point{d}}},
                 A::Vector{Vector{Point{d}}},
                 B::Vector{Vector{Point{d}}},
                 D::Vector{v},
                 J::Matrix{w}) where v<:AbstractVector{SMatrix{d,d,Float64,q}} where w<:AbstractVector{SMatrix{d,d,Float64,q}} where d where q

    n_m, n_M = size(J)
    N = length(F[1])
    for i=1:n_m
        for n=1:N
            f = F[i][n]
            x = X[i][n]
            a = A[i][n]
            b = B[i][n]
            z, z_f, z_x = PhiPhi(f,x,a,b)
            F[i][n] = z
            D[i][n] = z_f*D[i][n] + z_x
            for j=1:n_M
                J[i,j][n] = z_f*J[i,j][n]
            end
        end
    end
end

function PhiPhi!(F::Vector{<:AbstractVector{Point{d}}},
                 X::Vector{<:AbstractVector{Point{d}}},
                 A::Vector{<:AbstractVector{Point{d}}},
                 B::Vector{<:AbstractVector{Point{d}}},
                 J::Matrix{w}) where w<:AbstractVector{SMatrix{d,d,Float64,q}} where d where q

    n_m, n_M = size(J)
    N = length(F[1])
    for i=1:n_m
        for n=1:N
            f = F[i][n]
            x = X[i][n]
            a = A[i][n]
            b = B[i][n]
            z, z_f, z_x = PhiPhi(f,x,a,b)
            F[i][n] = z
            for j=1:n_M
                J[i,j][n] = z_f*J[i,j][n]
            end
        end
    end
end


function PhiPhi!(F::Vector{Vector{Point{d}}},
                X::Vector{Vector{Point{d}}},
                A::Vector{Vector{Point{d}}},
                B::Vector{Vector{Point{d}}},
                D::Vector{v}) where v<:AbstractVector{SMatrix{d,d,Float64,q}} where d where q

    n_m = size(D,1)
    N = length(F[1])
    for i=1:n_m
        for n=1:N
            f = F[i][n]
            x = X[i][n]
            a = A[i][n]
            b = B[i][n]
            z, z_f, z_x = PhiPhi(f,x,a,b)
            F[i][n] = z
            D[i][n] = z_f*D[i][n] + z_x
        end
    end
end


function PhiPhi(F::Vector{Vector{Point{d}}},X::Vector{Vector{Point{d}}},A::Vector{Vector{Point{d}}},B::Vector{Vector{Point{d}}},D::Vector{Vector{SMatrix{d,d,Float64,q}}}, J::Matrix{Vector{SMatrix{d,d,Float64,q}}}) where d where q
    FF = deepcopy(F)
    DD = deepcopy(D)
    JJ = deepcopy(J)
    PhiPhi!(FF,X,A,B,DD,JJ)
    return FF,DD,JJ
end
