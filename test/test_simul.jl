
horizon=40
seed=42


path = Pkg.dir("Dolo")

Pkg.build("QuantEcon")
import Dolo

filename = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model = Dolo.yaml_import(filename)
n_exp = 0
sigma = model.calibration.flat[:sig_z]

Dolo.simulation(model, sigma, n_exp, horizon, seed,zeros(0, 0))









# verbose=true
# verbose && @printf "%-8s%-10s%-10s%-10s%-5s\n" "t" model.symbols[:states][1] model.symbols[:states][2] model.symbols[:controls][1] model.symbols[:controls][2]
# verbose && println(repeat("-", 35))


#verbose && @printf "%-8s%-10s%-10s%-10s%-5s\n"  t round(s[1],2) round(s[2],2) round(x[1],2) round(x[2],2)
