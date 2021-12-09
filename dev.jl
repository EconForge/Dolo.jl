using Dolo 

module Temp
    using Dolo 

    struct DisG{ID, n_m, n_s, n_x, G_x, G_e}

        model::AbstractModel
        grid_endo
        grid_exo
        dprocess
        s0::Dolo.ListOfPoints{n_s}

        function DisG(model; 
            discretization=Dict())

            grid, dprocess = Dolo.discretize(model; discretization...)
            grid_endo = grid.endo
            grid_exo = grid.exo
        

            n_m = length(model.symbols[:exogenous])
            n_s = length(model.symbols[:states])
            n_x = length(model.symbols[:controls])

            s0 = Dolo.nodes(grid_endo)

            Gx = typeof(grid_exo)
            Ge = typeof(grid_endo)

            ID = id(model)

            new{ID, n_m, n_s, n_x, Gx, Ge}(model, grid_endo, grid_exo, dprocess, s0)
    
        end

    end
   
    function size_mu(disg::DisG)

        return length(Dolo.nodes(disg.grid_endo)) * length(Dolo.nodes(disg.grid_exo))
    
    end

    function ergodic_distribution(sol, disg::DisG)

        dr = sol.dr

        return Dolo.ergodic_distribution(disg.model, dr, disg.grid_exo, disg.grid_endo, disg.dprocess)[2]
    
    end


    function (G::DisG)(μ0,x0; exo = nothing)

        ∂G_∂μ = Dolo.new_transition(G.model, G.dprocess, x0, G.grid_exo, G.grid_endo; exo=exo)
        
        μ1 = Dolo.new_distribution(∂G_∂μ, μ0)

        ∂G_∂x = (x) -> μ0'* (ForwardDiff.gradient(new_transition, x))[3]

        return μ1, ∂G_∂μ, ∂G_∂x
    end

end

filename = "./examples/models/rbc_mc.yaml"
model = Model(filename)
sol = time_iteration(model)
dr = sol.dr
grid_endo, grid_exo = dr.grid_endo, dr.grid_exo
length(Dolo.nodes(grid_endo))
length(Dolo.nodes(grid_exo))
F = Dolo.Euler(model)
x0 = F.x0
μ0 = zeros(40)


G = Main.Temp.DisG(model)

Main.Temp.ergodic_distribution(sol, G)

Main.Temp.size_mu(G)


μ1, ∂G_∂μ, ∂G_∂x = G(μ0,x0)
∂G_∂μ
∂G_∂x

