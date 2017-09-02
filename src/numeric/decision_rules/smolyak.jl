@compat const SmolyakDR{S<:Grid,nx} = DecisionRule{S,<:SmolyakGrid,nx,Vector{Matrix{Float64}}}

#####
##### 1-argument decision rule
#####

function DecisionRule{nx}(grid_exo::EmptyGrid, grid_endo::SmolyakGrid, ::Union{Val{nx},Type{Val{nx}}})
    coeffs = [Array{Float64}(n_nodes(grid_endo), nx)]
    return DecisionRule(grid_exo, grid_endo, coeffs, Val{nx})
end

function set_values!(dr::SmolyakDR{<:EmptyGrid}, values::Vector{Matrix{Float64}})
    for i in 1:length(values)
        A_ldiv_B!(dr.itp[i], dr.grid_endo.B_nodes, values[i])
    end
end

function evaluate(dr::SmolyakDR{<:EmptyGrid}, z::AbstractMatrix)
    B = BM.evalbase(dr.grid_endo.smol_params, z)
    B*dr.itp[1]
end

####
#### UnstructuredGrid × CartesianGrid 2 continous arguments d.r.
####

## Smolyak!

function DecisionRule{nx}(grid_exo::UnstructuredGrid, grid_endo::SmolyakGrid, ::Union{Val{nx},Type{Val{nx}}})
    coeffs = [Array{Float64}(n_nodes(grid_endo), nx) for i in 1:n_nodes(grid_exo)]
    DecisionRule(grid_exo, grid_endo, coeffs, Val{nx})
end

function set_values!(dr::SmolyakDR{<:UnstructuredGrid}, values::Vector{Matrix{Float64}})
    for i in 1:length(values)
        A_ldiv_B!(dr.itp[i], dr.grid_endo.B_nodes, values[i])
    end
end

function evaluate(dr::SmolyakDR{<:UnstructuredGrid}, i::Int, z::AbstractMatrix)
    BM.evalbase(dr.grid_endo.smol_params, z) * dr.itp[i]
end
