using DataStructures

type DynareModelTestCase
    variables
    exogenous
    parameters
    equations
    calibration
    initval
    shocks
    modfile
end

# for fs2000.mod
fs_mod = let
    variables = ["m", "P", "c", "e", "W", "R", "k", "d", "n", "l", "gy_obs",
                 "gp_obs", "y", "dA"]
    exogenous = ["e_a", "e_m"]
    parameters = ["alp", "bet", "gam", "mst", "rho", "psi", "del"]
    equations = ["dA = exp(gam+e_a)",
                 "log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m",
                 "-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0",
                 "W = l/n",
                 "-(psi/(1-psi))*(c*P/(1-n))+l/n = 0",
                 "R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W",
                 "1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0",
                 "c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1)",
                 "P*c = m",
                 "m-1+d = l",
                 "e = exp(e_a)",
                 "y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a))",
                 "gy_obs = dA*y/y(-1)",
                 "gp_obs = (P/P(-1))*m(-1)/dA"]

    calibration = OrderedDict("alp" => "0.33",
                              "bet" => "0.99",
                              "gam" => "0.003",
                              "mst" => "1.011",
                              "rho" => "0.7",
                              "psi" => "0.787",
                              "del" => "0.02")
    initval = OrderedDict("k" => "6",
                          "m" => "mst",
                          "P" => "2.25",
                          "c" => "0.45",
                          "e" => "1",
                          "W" => "4",
                          "R" => "1.02",
                          "d" => "0.85",
                          "n" => "0.19",
                          "l" => "0.86",
                          "y" => "0.6",
                          "gy_obs" => "exp(gam)",
                          "gp_obs" => "exp(-gam)",
                          "dA" => "exp(gam)")

    # converted to vcov matrix -- you can choose a different convention if you want
    shocks = ["(0.014)^2" "0"
              "0"         "(0.005)^2"]

    modfile = Pkg.dir("Dolo", "examples", "models", "dynare", "fs2000.mod")
    DynareModelTestCase(variables, exogenous, parameters, equations,
                        calibration, initval, shocks, modfile)

end

rbc_mod = let
    variables = ["y", "c", "k", "i", "l", "y_l", "z"]
    exogenous = ["e"]
    parameters = ["beta", "psi", "delta", "alpha", "rho"]

    equations = ["(1/c) = beta*(1/c(+1))*(1+alpha*(k^(alpha-1))*(exp(z(+1))*l(+1))^(1-alpha)-delta)",
                 "psi*c/(1-l) = (1-alpha)*(k(-1)^alpha)*(exp(z)^(1-alpha))*(l^(-alpha))",
                 "c+i = y",
                 "y = (k(-1)^alpha)*(exp(z)*l)^(1-alpha)",
                 "i = k-(1-delta)*k(-1)",
                 "y_l = y/l",
                 "z = rho*z(-1)+e"]

    calibration = OrderedDict("alpha" => "0.33",
                              "beta" => "0.99",
                              "delta" => "0.023",
                              "psi" => "1.75",
                              "rho" => "0.95",
                              "sigma" => "(0.007/(1-alpha))")


    initval = OrderedDict("k" => "9",
                          "c" => "0.76",
                          "l" => "0.3",
                          "z" => "0",
                          "e" => "0")

    shocks = ["sigma^2"]
    modfile = Pkg.dir("Dolo", "examples", "models", "dynare", "rbc.mod")
    DynareModelTestCase(variables, exogenous, parameters, equations,
                        calibration, initval, shocks, modfile)
end
