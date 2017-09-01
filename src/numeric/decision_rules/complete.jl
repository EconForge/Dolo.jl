@compat const CompletePolyDR{S<:Grid,nx} = DecisionRule{S,<:RandomGrid,nx,Vector{Matrix{Float64}}}

#####
##### 1-argument decision rule
#####

function DecisionRule{n_s,nx}(grid_exo::EmptyGrid, grid_endo::RandomGrid{n_s}, ::Union{Val{nx},Type{Val{nx}}})
    coeffs = [Array{Float64}(BM.n_complete(n_s, 3), nx)]
    return DecisionRule(grid_exo, grid_endo, coeffs, Val{nx})
end

function set_values!(dr::CompletePolyDR{<:EmptyGrid}, values::Vector{Matrix{Float64}})
    B_grid = BM.complete_polynomial(nodes(dr.grid_endo), 3)
    for i in 1:length(values)
        A_ldiv_B!(dr.itp[i], B_grid, values[i])
    end
end

function evaluate(dr::CompletePolyDR{<:EmptyGrid}, z::AbstractMatrix)
    B = BM.complete_polynomial(z, 3)
    B*dr.itp[1]
end

####
#### UnstructuredGrid × CartesianGrid 2 continous arguments d.r.
####

function DecisionRule{n_s,nx}(grid_exo::UnstructuredGrid, grid_endo::RandomGrid{n_s}, ::Union{Val{nx},Type{Val{nx}}})
    coeffs = [
        Array{Float64}(BM.n_complete(n_s, 3), n_x)
        for i in 1:n_nodes(grid_exo)
    ]
    return DecisionRule(grid_exo, grid_endo, coeffs, Val{nx})
end

function set_values!(dr::CompletePolyDR{<:UnstructuredGrid}, values::Vector{Matrix{Float64}})
    B_grid = BM.complete_polynomial(nodes(dr.grid_endo), 3)
    for i in 1:length(values)
        A_ldiv_B!(dr.itp[i], B_grid, values[i])
    end
end

function evaluate(dr::CompletePolyDR{<:UnstructuredGrid}, i::Int, z::AbstractMatrix)
    B = BM.complete_polynomial(z, 3)
    B*dr.itp[i]
end
