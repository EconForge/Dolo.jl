import Dolo: get_grid, discretize
import Dolo: ListOfPoints, CachedDecisionRule

using Dolo
using StaticArrays

model = yaml_import("examples/models/rbc.yaml")

# @time sol = time_iteration(model;)



# @time soli = improved_time_iteration(model, sol.dr)




function discretize(model)
    grid = get_grid(model)
    dprocess = discretize(model.exogenous)
    return dprocess.grid, grid, dprocess
end

import Dolo: MSM

struct Euler{n_s, n_x, Gx, Ge}

    model::AbstractModel
    p:: SVector  ## SVector
    grid_endo
    grid_exo
    dprocess
    s0::ListOfPoints{n_s}
    x0::MSM{SVector{n_x, Float64}}
    dr::CachedDecisionRule
    cdr::CachedDecisionRule

    function Euler(model)

        p = SVector(model.calibration[:parameters]...)

        
        
        grid_exo, grid_endo, dprocess = discretize(model)
        
        x = SVector(model.calibration[:controls]...)
        
        N_m = Dolo.n_nodes(grid_exo)
        N_s = Dolo.n_nodes(grid_endo)

        n_s = length(model.symbols[:states])
        n_x = length(model.symbols[:controls])
        

        xxx = [[x for i=1:N_s] for j=1:N_m]

        x0 = MSM(xxx)

        s0 = Dolo.nodes(grid_endo)
        dr = Dolo.CachedDecisionRule(dprocess, grid_endo, xxx)
        cdr = Dolo.CachedDecisionRule(dprocess, grid_endo, xxx)


        Gx = typeof(grid_exo)
        Ge = typeof(grid_endo)

        new{n_s, n_x, Gx, Ge}(model, p, grid_endo, grid_exo, dprocess, s0, x0, dr, cdr)

    end

end


import Dolo: set_values!
import Dolo: euler_residuals

function (F::Euler)(x0::MSM, x1::MSM, set_future=true) 

    res = deepcopy(x0)

    if set_future
        set_values!(F.dr, x1)
    end
    
    rr =   euler_residuals(model,F.s0,x0,F.dr,F.dprocess,F.p;keep_J_S=false)

    return rr

end

# _,J_ij,S_ij =   euler_residuals(model,s,x,ddr,dprocess,p,keep_J_S=true,set_dr=true)



@time F = Euler(model);



z0 = F.x0;
F(z0,z0);

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


rrrr = F(z0,z0,true);



# import ProfileView
# ProfileView.@profview res = F(z0,z0,true);




@time res = F(z0,z0,true);


@time res = F(z0,z0,false);


function df_A(F, z0, z1; set_future=false)

    fun  = z->F(z, z1, false)

    J = Dolo.DiffFun(fun, z0)[2]

    return (J)

end

import Dolo: LinearThing

function df_B(F, z0, z1; set_future=false)

    ddr_filt = F.cdr

    # res = deepcopy(z0)

    if set_future
        set_values!(F.dr, z1)
    end

    _,J_ij,S_ij  =   euler_residuals(model, F.s0, z0 , F.dr, F.dprocess, F.p; keep_J_S=true)

    L = LinearThing(J_ij, S_ij, ddr_filt)

    return L

end


J = df_A(F, z0, z0);

L = df_B(F, z0, z0);

function *(x::MSM, y::MSM)
    data = x.data .* y.data
    return MSM(data, x.sizes)
end

# function *(L::LinearThing,x::MSM{SVector{n_x, Float64}}) where n_x
#     res = L*x.views
#     return MSM(cat(res...; dims=1), x.sizes)
# end

function prediv!(L::LinearThing,x::MSM)
    for i=1:size(L.M_ij,1)
        N = length(x.views[i])
        for j=1:size(L.M_ij,2)
            for n=1:N
                L.M_ij[i,j][n] = x.views[i][n] \ L.M_ij[i,j][n]
            end
        end

    end
    return L
end

function premult!(L::LinearThing,x::MSM)
    for i=1:size(L.M_ij,1)
        N = length(x.views[i])
        for j=1:size(L.M_ij,2)
            for n=1:N
                L.M_ij[i,j][n] = x.views[i][n] * L.M_ij[i,j][n]
            end
        end

    end
    return L
end

function mult!(L::LinearThing,x::Number)
    for i=1:size(L.M_ij,1)
        for j=1:size(L.M_ij,2)
            L.M_ij[i,j][:] *= x
        end

    end
end

function invert!(x::MSM)
    for i=1:length(x.data)
        x.data[i] = inv(x.data[i] )
    end
end



premult!(L, J)

function ldiv(a::MSM{S}, b::MSM{T}) where S where T
    r = a.data .\ b.data
    MSM(r, a.sizes)
end

# @code_warntype df_A(F, z0, z0);

import Dolo: newton


function loop_ti(F::Euler, z0::MSM{V}; tol_η=1e-7, tol_ε=1e-8, T=100) where V

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


        δ = z1 - z0

        err = norm(δ)

        if err<tol_η
            return 
        end
        z0 = z1

    end

    return err
end

@time loop_ti(F, z0);


function loop_iti(F::Euler, z0::MSM{V}; verbose=true, tol_η=1e-7, tol_ε=1e-8, tol_κ=1e-8, T=500, K=500, switch=5, mode=:iti) where V

    local err

    for t=1:T

        r = F(z0, z0, true)

        err_1 = norm(r)


        if mode==:newton

            sol = newton(u-> F(u,z0,false), z0)
            z1 = sol.solution
    
            δ = z1 - z0
            count = sol.iterations
    
        else

            J = df_A(F, z0, z0 ; set_future=false)

            if t <=switch

                # incomplete time iteration step
                δ = -J\r

                count = 0

            else

                L = df_B(F, z0, z0 ; set_future=false)

                mult!(L, -1.0)
                prediv!(L, J) # J\L

                π = -J\r

                count = 0

                u = π
                δ = π
                for i=1:K
                    count +=1
                    u = L*u
                    δ += u
                    if norm(u)<tol_κ
                        break
                    end
                end
            
            end

        end

        err = norm(δ)

        if verbose
            println("$t | $err_1 | $err | $count")
        end

        if err<tol_η
            return 
        end
        z0 = z0 + δ


    end

    return err
end

@time res =  loop_iti(F, z0; T=500, switch=0, K=50);

@time res =  loop_iti(F, z0; T=1000, switch=500, K=50);

@time res =  loop_iti(F, z0; T=1000, verbose=false, switch=500, K=50, mode=:newton);

@time time_iteration(model, verbose=false, m)

power_iteration(L)



function power_iteration(L)
    N = length(L.M_ij)*length(L.M_ij[1])
    z0 = rand(500)

    z0 = z0./maximum(abs, z0)
    for k = 1:100
        z1 = L*z0
        n1 = maximum(abs, z1)
        z0 = z1/n1
        println("norm: $(n1)")
    end
end



@time sol = time_iteration(model; verbose=false, complementarities=false);



using ProfileView

@time res =  loop_iti(F, z0; T=500, switch=0, K=50, verbose=false);
@profview  loop_iti(F, z0; T=500, switch=0, K=50, verbose=false);





N = 1000000

m0, s0, x0, p = [SVector(e...) for e in model.calibration[:exogenous,:states,:controls,:parameters]]


mvv = [m0 for i in 1:N]
svv = [s0 for i in 1:N]
xvv = [x0 for i in 1:N]
pvv = [p for i in 1:N]


out = copy(xvv)
out[:]*=0

@time arbitrage(model, m0, s0, x0, m0, s0, x0, p);

@time arbitrage( model, mvv, svv, xvv, mvv, svv, xvv, pvv, out);


fun(u...) = arbitrage(model, u...);

fun(m0, s0, x0,m0, s0, x0, p)

@time fun(mvv, svv, xvv, mvv, svv, xvv, pvv);

@time fun.(mvv, svv, xvv, mvv, svv, xvv, pvv);

fff = model.factories[:transition]

import Dolang
code = Dolang.gen_kernel(fff, [0]; funname=:transit)

transit = eval(code)


@time Dolo.transition(model, mvv, svv, xvv, mvv, pvv);

@time transit.(mvv, svv, xvv, mvv, pvv);


function transit(m::SVector{1, Float64}, s::SVector{1, Float64}, x::SVector{2, Float64}, M::SVector{1, Float64}, p::SVector{9, Float64})
    z_m1_ = m[1]
    k_m1_ = s[1]
    n_m1_ = x[1]
    i_m1_ = x[2]
    z__0_ = M[1]
    beta_ = p[1]
    sigma_ = p[2]
    eta_ = p[3]
    chi_ = p[4]
    delta_ = p[5]
    alpha_ = p[6]
    rho_ = p[7]
    zbar_ = p[8]
    sig_z_ = p[9]
    y_m1_ = (exp(z_m1_) * k_m1_ ^ alpha_) * n_m1_ ^ (1.0 - alpha_)
    w_m1_ = ((1.0 - alpha_) * y_m1_) / n_m1_
    rk_m1_ = (alpha_ * y_m1_) / k_m1_
    k__0_ = (1.0 - delta_) * k_m1_ + i_m1_
    c_m1_ = y_m1_ - i_m1_
    oo_1 = SVector(k__0_)
    res_ = (oo_1,)
    return res_
end

function transit2(m::SVector{1, Float64}, s::SVector{1, Float64}, x::SVector{2, Float64}, M::SVector{1, Float64}, p::SVector{9, Float64})
    # z_m1_ = m[1]
    k_m1_ = s[1]
    # n_m1_ = x[1]
    i_m1_ = x[2]
    # z__0_ = M[1]
    # beta_ = p[1]
    # sigma_ = p[2]
    # eta_ = p[3]
    # chi_ = p[4]
    delta_ = p[5]
    # alpha_ = p[6]
    # rho_ = p[7]
    # zbar_ = p[8]
    # sig_z_ = p[9]
    # y_m1_ = (exp(z_m1_) * k_m1_ ^ alpha_) * n_m1_ ^ (1.0 - alpha_)
    # w_m1_ = ((1.0 - alpha_) * y_m1_) / n_m1_
    # rk_m1_ = (alpha_ * y_m1_) / k_m1_
    k__0_ = (1.0 - delta_) * k_m1_ + i_m1_
    # c_m1_ = y_m1_ - i_m1_
    oo_1 = SVector(k__0_)
    res_ = (oo_1,)
    return res_
end




function transit_stupid(z_m1_, k_m1_, n_m1_, i_m1_, z__0_, beta_, sigma_, eta_, chi_, delta_, alpha_, rho_, zbar_, sig_z_)
    begin
        
    y_m1_ = (exp(z_m1_) * k_m1_ ^ alpha_) * n_m1_ ^ (1.0 - alpha_)
    w_m1_ = ((1.0 - alpha_) * y_m1_) / n_m1_
    rk_m1_ = (alpha_ * y_m1_) / k_m1_
    k__0_ = (1.0 - delta_) * k_m1_ + i_m1_

    end
    
    c_m1_ = y_m1_ - i_m1_
    oo_1 = SVector(k__0_)
    res_ = (oo_1,)
    return res_

end

function transit_stupid(m::SVector{1, Float64}, s::SVector{1, Float64}, x::SVector{2, Float64}, M::SVector{1, Float64}, p::SVector{9, Float64})
    transit_stupid(m..., s..., x..., M..., p...)
end

transit_stupid(m0, s0, x0, m0, p)

@time transit_stupid.(mvv, svv, xvv, mvv, pvv);


@time transit.(mvv, svv, xvv, mvv, pvv);

@time transit2.(mvv, svv, xvv, mvv, pvv);


@code_warntype transit(m0, s0, x0, m0, p);
