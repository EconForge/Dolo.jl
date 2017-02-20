##function QuantEcon.simulate(m::AbstractNumericModel, dr::AbstractDecisionRule,
#                            s0::AbstractVector=m.calibration[:states];
#                            n_exp::Int=0, horizon::Int=40,
#                            seed::Int=42,
#                            forcing_shocks::AbstractMatrix=zeros(0, 0))
#
#
#forcing_shocks=zeros(0, 0)
#verbose::Bool=true
horizon=40
seed=42
#irf = n_exp == 0 ? true : false

path = Pkg.dir("Dolo")

Pkg.build("QuantEcon")
import Dolo

filename = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model = Dolo.yaml_import(filename)
n_exp = 0
n_exp = max(n_exp, 1)

# extract data from model
calib = model.calibration
params = calib[:parameters]
#ny = length(calib[:auxiliaries])
#has_aux = ny > 0

sigma = model.calibration.flat[:sigma]

# calculate initial controls using decision rule
s0=model.calibration[:states]
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)
dr
x0 = dr(s0)

# get number of states and controls
ns = length(s0)
nx = length(x0)
nsx = nx+ns

# TODO: talk to Pablo and make a decision about this
s_simul = Array(Float64, n_exp, ns, horizon)
x_simul = Array(Float64, n_exp, nx, horizon)
for i in 1:n_exp
  s_simul[i, :, 1] = s0
  x_simul[i, :, 1] = x0
end
# NOTE: this will be empty if ny is zero. That's ok. Our call to `cat`  #       below will work either way  y_simul = Array(Float64, n_exp, ny, horizon)


srand(seed)

using Distributions
d = Normal(0,sigma)
epsilons = rand(d,horizon)
# not working
#ϵ_dist = Dolo.MvNormal(sigma)
#epsilons = rand(ϵ_dist, n_exp)'

#verbose=true
#verbose && @printf "%-8s%-10s%-10s%-10s%-5s\n" "t" model.symbols[:states][1] model.symbols[:states][2] model.symbols[:controls][1] model.symbols[:controls][2]
#verbose && println(repeat("-", 35))


for t in 2:horizon
    #if irf
    #  if !isempty(forcing_shocks) && t < size(forcing_shocks, 2)
    #      epsilons = forcing_shocks[t, :]'
    #  else
    #      epsilons = zeros(size(sigma, 1), 1)
    #  end
    ##else
    ##    epsilons = rand(ϵ_dist, n_exp)'
    #end

    #s = view(s_simul, :, :, t)
    s = copy(view(s_simul, :, :, t))
    x = copy(view(x_simul, :, :, t))
    x = dr(s)
    x_simul[:, :, t] = x

    if t < horizon
      ss = copy(view(s_simul, :, :, t+1))

      ss = Dolo.transition!(model, s[1:ns],  [epsilons[1]], s[1:ns], x[1:nx], [epsilons[t+1]], params)

      s_simul[:, :, t+1] = ss
      # verbose && @printf "%-8s%-10s%-10s%-10s%-5s\n"  t round(ss[1],2) round(ss[2],2) round(xx[1],2) round(xx[2],2)
    end


end

cat(2, s_simul, x_simul)::Array{Float64,3}
