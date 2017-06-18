function filter_coeffs(a::Vector{Float64}, b::Vector{Float64}, orders::Vector{Int64}, data::Array{Float64})

    dinv = (b - a)./(orders.-1)
    d = length(a)
    if d == 1
        coeffs = filter_coeffs_1d(dinv,data)
    elseif d == 2
        coeffs = filter_coeffs_2d(dinv,data)
    elseif d == 3
        coeffs = filter_coeffs_3d(dinv,data)
    elseif d == 4
        coeffs = filter_coeffs_4d(dinv,data)
    else
        error()
    end
    return coeffs
end

function ff(a::Vector{Float64}, b::Vector{Float64}, orders::Vector{Int64}, data::Array{Float64})

    dinv = (b - a)./(orders.-1)
    d = length(a)
    if d == 1
        coeffs = filter_coeffs_1d(dinv,data)
    elseif d == 2
        coeffs = filter_coeffs_2d(dinv,data)
    elseif d == 3
        coeffs = filter_coeffs_3d(dinv,data)
    elseif d == 4
        coeffs = filter_coeffs_4d(dinv,data)
    else
        error()
    end
    return coeffs
end

# function find_coeffs_1d_old(delta_inv, M, data)
#
#     # args: double, int, double[:], double[:
#
#     N = M + 2
#
#     data = data[:]
#
#     if length(data) == N
#         rhs = data
#     else
#         rhs = [0; data; 0]
#     end
#
#     basis = [1.0/6.0, 2.0/3.0, 1.0/6.0]
#
#     vals = repeat( basis, outer=[M])
#     co_x = repeat( collect(2:M+1), inner=[3])
#     co_y = repeat( collect(1:3) , outer=[M]) + co_x .- 2
#
#     db = 4
#     initial = [1.0,-2.0,1.0,0.0]*delta_inv*delta_inv
#     final = [0.0,1.0,-2.0,1.0]*delta_inv*delta_inv
#
#     vals = [initial; vals; final]
#     co_x = [ones(Int,db); co_x; ones(Int,db)*(M+2)]
#     co_y = [1:db; co_y; M+3-db:M+2]
#
#     spmat = sparse(co_x, co_y, vals)
#
#     sol = spmat \ rhs
#
#     return sol
#
# end

function solve_deriv_interp_1d!(bands::Matrix{Float64}, coefs::Vector{Float64})

    M = size(coefs,1) - 2

    # Solve interpolating equations
    # First and last rows are different

    bands[1, 2] /= bands[1, 1]
    bands[1, 3] /= bands[1, 1]
    bands[1, 4] /= bands[1, 1]
    bands[1, 1] = 1.0
    bands[2, 2] -= bands[2, 1] * bands[1, 2]
    bands[2, 3] -= bands[2, 1] * bands[1, 3]
    bands[2, 4] -= bands[2, 1] * bands[1, 4]
    bands[1, 1] = 0.0
    bands[2, 3] /= bands[2, 2]
    bands[2, 4] /= bands[2, 2]
    bands[2, 2] = 1.0

    # Now do rows 2 through M+1
    for row = 3:(M+1)
        bands[row, 2] -= bands[row, 1] * bands[row - 1, 3]
        bands[row, 4] -= bands[row, 1] * bands[row - 1, 4]
        bands[row, 3] /= bands[row, 2]
        bands[row, 4] /= bands[row, 2]
        bands[row, 1] = 0.0
        bands[row, 2] = 1.0
    end
    # Do last row
    bands[M + 2, 2] -= bands[M + 2, 1] * bands[M, 3]
    bands[M + 2, 4] -= bands[M + 2, 1] * bands[M, 4]
    bands[M + 2, 3] -= bands[M + 2, 2] * bands[M + 1, 3]
    bands[M + 2, 4] -= bands[M + 2, 2] * bands[M + 1, 4]
    bands[M + 2, 4] /= bands[M + 2, 3]
    bands[M + 2, 3] = 1.0

    coefs[M + 2] = bands[M + 2, 4]
    # Now back substitute up
    # for row in (M+1):2
    for ii in 1:M
        row = M+2 - ii
        coefs[row] = bands[row, 4] - bands[row, 3] * coefs[row + 1]
    end
    # Finish with first row
    coefs[1] = bands[1, 4] - bands[1, 2] * coefs[2] - bands[1, 3] * coefs[3]

end


function find_coeffs_1d!(delta_inv, M, data, coefs)


    basis = [1.0/6.0, 2.0/3.0, 1.0/6.0]
    bands = zeros(M+2, 4)

    # Setup boundary conditions
    abcd_left = zeros(4)
    abcd_right = zeros(4)

    # Left boundary
    abcd_left[1] = 1.0 * delta_inv * delta_inv
    abcd_left[2] = -2.0 * delta_inv * delta_inv
    abcd_left[3] = 1.0 * delta_inv * delta_inv
    abcd_left[4] = 0

    # Right boundary
    abcd_right[1] = 1.0 * delta_inv * delta_inv
    abcd_right[2] = -2.0 * delta_inv * delta_inv
    abcd_right[3] = 1.0 * delta_inv * delta_inv
    abcd_right[4] = 0

    for i in 1:4
        bands[1, i] = abcd_left[i]
        bands[M + 2, i] = abcd_right[i]
    end
    for i in 2:(M+1)
        for j in 1:3
            bands[i, j] = basis[j]
            bands[i, 4] = data[i-1]
        end
    end

    solve_deriv_interp_1d!(bands, coefs)

    return coefs
end

function find_coeffs_1d(delta_inv, M, data)
    coefs = zeros(M+2)
    find_coeffs_1d!(delta_inv, M, data, coefs)
    return coefs
end

function filter_coeffs_1d(dinv, data)

  M = size(data,1)
  N = M+2

  coefs = find_coeffs_1d(dinv[1], M, data)

  return coefs

end

function filter_coeffs_2d(dinv, data)

    Mx = size(data,1)
    My = size(data,2)

    Nx = Mx+2
    Ny = My+2

    coefs = zeros(Nx,Ny)


    # First, solve in the X-direction
    for iy in 1:My
        coefs[:,iy+1] = find_coeffs_1d(dinv[1], Mx, data[:,iy])
    end

    # Now, solve in the Y-direction
    for ix in 1:Nx
        coefs[ix,:] = find_coeffs_1d(dinv[2], My, coefs[ix,:])
    end

    return coefs

end

function filter_coeffs_3d(dinv, data)

    Mx = size(data,1)
    My = size(data,2)
    Mz = size(data,3)

    Nx = Mx+2
    Ny = My+2
    Nz = Mz+2

    coefs = zeros(Nx,Ny,Nz)


    # First, solve in the X-direction
    for iy in 1:My
        for iz in 1:Mz
            coefs[:,iy+1,iz+1] = find_coeffs_1d(dinv[1], Mx, data[:,iy,iz])
        end
    end

    # Now, solve in the Y-direction
    for ix in 1:Nx
        for iz in 1:Mz
            coefs[ix,:,iz+1] = find_coeffs_1d(dinv[2], My, coefs[ix,:,iz+1])
        end
    end

    # Now, solve in the Z-direction
    for ix in 1:Nx
        for iy in 1:Ny
            coefs[ix,iy,:] = find_coeffs_1d(dinv[3], Mz, coefs[ix,iy,:])
        end
    end

    return coefs

end

function filter_coeffs_4d(dinv, data)

    Mx = size(data,1)
    My = size(data,2)
    Mz = size(data,3)
    Mz4 = size(data,4)

    Nx = Mx+2
    Ny = My+2
    Nz = Mz+2
    Nz4 = Mz4+2

    coefs = zeros(Nx,Ny,Nz,Nz4)


    # First, solve in the X-direction
    for iy in 1:My
        for iz in 1:Mz
            for iz4 in 1:Mz4
                coefs[:,iy+1,iz+1,iz4+1] = find_coeffs_1d(dinv[1], Mx, data[:,iy,iz,iz4])
            end
        end
    end

    # Now, solve in the Y-direction
    for ix in 1:Nx
        for iz in 1:Mz
            for iz4 in 1:Mz4
                coefs[ix,:,iz+1,iz4+1] = find_coeffs_1d(dinv[2], My, coefs[ix,:,iz+1,iz4+1])
            end
        end
    end

    # Now, solve in the Z-direction
    for ix in 1:Nx
        for iy in 1:Ny
            for iz4 in 1:Mz4
                coefs[ix,iy,:,iz4+1] = find_coeffs_1d(dinv[3], Mz, coefs[ix,iy,:,iz4+1])
            end
        end
    end
    # Now, solve in the Z4-direction
    for ix in 1:Nx
        for iy in 1:Ny
            for iz in 1:Mz
                coefs[ix,iy,iz,:] = find_coeffs_1d(dinv[4], Mz, coefs[ix,iy,iz,:])
            end
        end
    end

    return coefs

end
