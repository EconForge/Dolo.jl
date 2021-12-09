using Dolo 

module Temp
    using Dolo 

    struct disG{ID, n_m, n_s, n_x, G_x, G_e}

        model::AbstractModel
        grid_endo
        grid_exo
        dprocess
        s0::Dolo.ListOfPoints{n_s}

        function disG(model; 
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
   
    function size_mu(DISG::disG)

        return length(Dolo.nodes(DISG.grid_endo)) * length(Dolo.nodes(DISG.grid_exo))
    
    end

    function ergodic_distribution(DISG::disG)

        sol = time_iteration(DISG.model)
        dr = sol.dr

        return Dolo.ergodic_distribution(DISG.model, dr, DISG.grid_exo, DISG.grid_endo, DISG.dprocess)[2]
    
    end

    function G(μ0,x0, DISG::disG; exo = nothing)

        sol = time_iteration(DISG.model)

        ∂G_∂μ = Dolo.new_transition(DISG.model, sol, x0, DISG.grid_exo, DISG.grid_endo; exo=exo)
        
        μ1 = Dolo.new_distribution(∂G_∂μ, μ0)

        ∂G_∂x = (x) -> (ForwardDiff.gradient(new_transition, x).*μ0)[3]

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


disg_struct = Main.Temp.disG(model)

Main.Temp.ergodic_distribution(disg_struct)

Main.Temp.size_mu(disg_struct)


μ1, ∂G_∂μ, ∂G_∂x = Main.Temp.G(μ0,x0, disg_struct)
∂G_∂μ
∂G_∂x
