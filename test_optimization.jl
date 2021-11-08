using StaticArrays


N = 100000000

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

    return k__0_
    # oo_1 = SVector(k__0_)
    # res_ = (oo_1,)
    # return res_
end

function transit2(m::SVector{1, Float64}, s::SVector{1, Float64}, x::SVector{2, Float64}, M::SVector{1, Float64}, p::SVector{9, Float64})                                                                                                                             
    k_m1_ = s[1]                                                                                                                                                                                                                                                                   
    i_m1_ = x[2]                                                                                                                                                                                                                                                                   

    delta_ = p[5]
    k__0_ = (1.0 - delta_) * k_m1_ + i_m1_
    return k__0_
    # oo_1 = SVector(k__0_)
    # res_ = (oo_1,)
    return res_
end



m0 = SVector(0.0)
s0 = SVector(9.5)
x0 = SVector(0.33, 0.23)
p = SVector(0.96, 2.00, 0.4, 0.5,  0.1, 0.3, 0.9, 0.0, 0.01)

transit(m0, s0, x0, m0, p)

mv = [m0 for i=1:N]
sv = [s0 for i=1:N]
xv = [x0 for i=1:N]
pv = [p  for i=1:N]


println("Warmup")
@time transit.(mv, sv, xv, mv, pv);
@time transit2.(mv, sv, xv, mv, pv);


println("Real Test")

@time transit.(mv, sv, xv, mv, pv);
@time transit2.(mv, sv, xv, mv, pv);
