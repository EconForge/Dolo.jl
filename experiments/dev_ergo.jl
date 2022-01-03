using Dolo


model = yaml_import("examples/models/rbc_mc.yaml")

sol = Dolo.improved_time_iteration(model)

Π = Dolo.transition_matrix(model, sol)
dΠ = Dolo.transition_matrix(model, sol; diff=true)



module Temp
    
    using Dolo: AbstractModel, ListOfPoints, MSM, n_nodes
    using Dolo: id, transition_matrix
    using StaticArrays
    struct distG{ID, n_m, n_s, n_x, Gx, Ge}

        model::AbstractModel

        grid_endo
        grid_exo
        dprocess

        s0::ListOfPoints{n_s}
        x0::MSM{SVector{n_x, Float64}}

    end

    function distG(model, sol)
        ID = id(model)
        grid_endo = sol.dr.grid_endo
        grid_exo = sol.dr.grid_exo
        dprocess = sol.dprocess
        s0 = grid_endo.nodes
        x0 = MSM([sol.dr(i, s0) for i=1:max(1,n_nodes(grid_exo))])
        n_m = ndims(grid_exo)
        n_s = ndims(grid_endo)
        n_x = length(x0.data[1])
        Gx = typeof(grid_exo)
        Ge = typeof(grid_endo)
        return distG{ID, n_m, n_s, n_x, Gx, Ge}(model, grid_endo, grid_exo, dprocess, s0, x0)
    end

    function (G::distG)(mu0, x0)

        P = transition_matrix(G.model, G.dprocess, x0, G.grid_exo, G.grid_endo; exo=nothing, diff=false)
        mu1 = P'*mu0
    end

end

G = Temp.distG(model, sol)

mu0 = ones( Dolo.n_nodes(G.grid_endo)*Dolo.n_nodes(G.grid_exo))
x0 = G.x0

@time P = G(mu0, G.x0);
