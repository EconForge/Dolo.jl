import Dolo: get_grid, discretize
import Dolo: ListOfPoints, CachedDecisionRule

using Dolo
using StaticArrays

model = yaml_import("examples/models/rbc.yaml")

model = yaml_import("examples/models/rbc_iid.yaml")

mvn = Dolo.MvNormal([0.2 0.1; 0.1 0.2])

𝔼 = Dolo.Euler(model)

𝔼.

@time sol = time_iteration(model;)

@time res =  Dolo.loop_ti(model; T=500);

@time res =  Dolo.loop_iti(model; T=500, switch=0, K=30);

@time res =  Dolo.loop_iti(model; T=1000, switch=5, K=50);

@time res =  loop_iti(F, z0; T=1000, verbose=false, switch=500, K=50, mode=:newton);

# @time soli = improved_time_iteration(model, sol.dr)




@time time_iteration(model, verbose=false)

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





N = 100000

m0, s0, x0, p = [SVector(e...) for e in model.calibration[:exogenous,:states,:controls,:parameters]]


mvv = [m0 for i in 1:N]
svv = [s0 for i in 1:N]
xvv = [x0 for i in 1:N]
pvv = [p for i in 1:N]


out = copy(xvv)
out[:]*=0

@time arbitrage(model, m0, s0, x0, m0, s0, x0, p);

@time arbitrage( model, mvv, svv, xvv, mvv, svv, xvv, pvv);

fun(u...) = arbitrage(model, u...);

@time @tturbo fun.( mvv, svv, xvv, mvv, svv, xvv, pvv);


@code_warntype arbitrage(model, m0, s0, x0, m0, s0, x0, p);

@using 

fun(m0, s0, x0,m0, s0, x0, p)

@time fun(mvv, svv, xvv, mvv, svv, xvv, pvv);

@time fun.(mvv, svv, xvv, mvv, svv, xvv, pvv);

fff = model.factories[:arbitrage]

import Dolang
code = Dolang.gen_kernel(fff, [0]; funname=:arbitrage)

transit = eval(code)


@time Dolo.transition(model, mvv, svv, xvv, mvv, pvv);

@time transit.(mvv, svv, xvv, mvv, pvv);


@fastmath function transit0(m::SVector{1, Float64}, s::SVector{1, Float64}, x::SVector{2, Float64}, M::SVector{1, Float64}, p::SVector{9, Float64})
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
    k__0_ = (1.0 - delta_) * k_m1_ ^ i_m1_
    c_m1_ = y_m1_ - i_m1_
    oo_1 = SVector(k__0_)
    res_ = (oo_1,)
    return res_
end

@fastmath function transit2(m::SVector{1, Float64}, s::SVector{1, Float64}, x::SVector{2, Float64}, M::SVector{1, Float64}, p::SVector{9, Float64})
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

    k__0_ = (1.0 - delta_) * k_m1_ ^i_m1_
    oo_1 = SVector(k__0_)
    res_ = (oo_1,)
    return res_
end



@code_llvm transit0(m0, s0, x0, m0, p)

@code_llvm transit2(m0, s0, x0, m0, p)




@time transit0.(mvv, svv, xvv, mvv, pvv);

@time transit2.(mvv, svv, xvv, mvv, pvv);


@time  transition.(mvv, svv, xvv, mvv, pvv);

tchuk

@code_warntype transit(m0, s0, x0, m0, p);




fff = model.factories[:transition]



@code_llvm sin(0.1)

function unaccessed_variable_1(m,s)

    m_1 = m[1]
    s_1 = s[1]
    x = m_1 + s_1 # s_2 is unused

end



function unaccessed_variable_2(m,s)

    m_1 = m[1]
    s_1 = s[1]
    s_2 = s[2]
    x = m_1 + s_1 # s_2 is unused
    return x

end


function unaccessed_variable_3(m,s)

    m_1 = m[1]
    s_1 = s[1]
    s_2 = s[2]
    z = s_2 ^ m_1 # z is unused
    x = m_1 + s_1 # s_2 is unused
    return x

end

function unaccessed_variable_3(m, s)

    m_1 = m[1]
    s_1 = s[1]
    s_2 = s[2]
    z = s_2 ^ m_1 # z is unused
    x = m_1 + s_1 # s_2 is unused
    return x
end
function unaccessed_variable_3(m::SVector{1,Float64},s::SVector{2,Float64})

    m_1 = m[1]
    s_1 = s[1]
    s_2 = s[2]
    z = s_2 ^ m_1 # z is unused
    x = m_1 + s_1 # s_2 is unused
    return x
end

a = [0.0]
b = [0.1, 0.2]

@code_llvm @fastmath unaccessed_variable_1(a, b)
@code_llvm @fastmath unaccessed_variable_2(a, b)
@code_llvm @fastmath unaccessed_variable_3(a, b)



av = SVector(0.0)
bv = SVector(0.1, 0.2)

@code_llvm @fastmath unaccessed_variable_1(a, b)
@code_llvm @fastmath unaccessed_variable_2(a, b)
@code_llvm @fastmath unaccessed_variable_3(a, b)


@code_llvm @fastmath unaccessed_variable_3(av, bv)



### import Dolang

import Dolang

code = Dolang.gen_kernel(model.factories[:transition], [0,1])

using MacroTools

clean_unused(code) 
eval( clean_unused(code) )