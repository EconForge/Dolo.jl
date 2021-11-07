import Dolo: get_grid, discretize

using Dolo
using StaticArrays

model = yaml_import("examples/models/rbc.yaml")

@time sol = time_iteration(model;)



@time soli = improved_time_iteration(model, sol.dr)




function discretize(model)
    grid = get_grid(model)
    dprocess = discretize(model.exogenous)
    return dprocess.grid, grid, dprocess
end

import Dolo: MSM

struct Euler

    model::AbstractModel
    p  ## SVector
    grid_endo
    grid_exo
    dprocess
    s0
    x0 ::MSM
    dr

    function Euler(model)

        p = SVector(model.calibration[:parameters]...)

        
        
        grid_exo, grid_endo, dprocess = discretize(model)
        
        x = SVector(model.calibration[:controls]...)
        
        N_m = Dolo.n_nodes(grid_exo)
        N_s = Dolo.n_nodes(grid_endo)
        n_x = length(model.symbols[:controls])
        

        xxx = [[x for i=1:N_s] for j=1:N_m]

        x0 = MSM(xxx)

        s0 = Dolo.nodes(grid_endo)
        dr = Dolo.CachedDecisionRule(dprocess, grid_endo, xxx)
        

        new(model, p, grid_endo, grid_exo, dprocess, s0, x0, dr)

    end

end


import Dolo: set_values!

function (F::Euler)(x0::MSM, x1::MSM, set_future=true) 

    res = deepcopy(x0)

    if set_future
        set_values!(F.dr, x1)
    end
    
    rr = Dolo.euler_residuals_ti!(res, F.model, F.dprocess, F.s0, x0, F.p, F.dr) 

    return (rr)

end


@time F = Euler(model);

z0 = F.x0;

import Base: -, \, +, /, *

function -(a::MSM)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = -a.data
    return MSM(data, sizes)
end


function -(a::MSM, b::MSM)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data .- b.data
    return MSM(data, sizes)
end



function +(a::MSM, b::MSM)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data .+ b.data
    return MSM(data, sizes)
end


function \(a::MSM, b::MSM)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data .\ b.data
    return MSM(data, sizes)
end

function /(a::MSM, b::MSM)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data ./ b.data
    return MSM(data, sizes)
end

function /(a::MSM, b::Number)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data ./ b
    return MSM(data, sizes)
end

function *(a::MSM, b::Number)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data .* b
    return MSM(data, sizes)
end



import Dolo: maxabs

norm(a::MSM) = maximum( u-> maximum(abs,u), a.data )
maxabs(a::MSM) = maximum( u-> maximum(abs,u), a.data )


F(z0,z0,true);


@time res = F(z0,z0,true);


@time res = F(z0,z0,false);

import Dolo: vecvec

function df_A(F, z0, z1; set_future=false)

    fun  = z->vecvec(F(MSM(z), z1, false))

    zz0 = vecvec(z0)

    J = Dolo.DiffFun(fun, zz0)[2]

    return MSM(J)

end
# df_B(F, x0, x1)


@time J = df_A(F, z0, z0);

import Dolo: newton


function loop(F::Euler, z0::MSM{V}; tol_η=1e-7, tol_ε=1e-8, T=100) where V

    z1 = deepcopy(z0)

    local err
    for t=1:T

        res = F(z0, z0, true)
        err_1 = norm(res)
        if err_1<tol_ε
            return
        end


        sol = newton(u-> F(u,z0,false), z0)

        z1 = sol.solution
        # J = df_A(F, z0, z0)
        # δ = - (J\res)

        δ = z1 - z0

        err = norm(δ)

        if err<tol_η
            return 
        end
        z0 = z1
        # println(err)

    end

    return err
end

@time loop(F, z0);






@time time_iteration(model; verbose=false, complementarities=false);

