using Dolo 
using BenchmarkTools
using StaticArrays
using FiniteDiff


model = yaml_import("examples/models/rbc.yaml")
sol = Dolo.improved_time_iteration(model)


G = Dolo.distG(model, sol)
z10 = SVector(model.calibration[:exogenous]...) 
z20 = z10 
x0 = G.x0
x0_flat = cat(G.x0.data...; dims=1)
μ0 = G.μ0



μ1, ∂G_∂μ, ∂G_∂x, ∂G_∂z1, ∂G_∂z2 = G(μ0, x0; exo = [z10,z20], diff = true)

Jμ_exact = convert(Matrix, ∂G_∂μ)
Jμ_num = FiniteDiff.finite_difference_jacobian(mu -> G(mu, x0, exo = [z10,z20]), μ0)

Jx_exact = convert(Matrix, ∂G_∂x)
Jx_num = FiniteDiff.finite_difference_jacobian(x -> G(μ0, x; exo = [z10,z20]), x0_flat)

Jz1_exact = convert(Matrix, ∂G_∂z1)
Jz1_num = FiniteDiff.finite_difference_jacobian(z1 -> G(μ0, x0; exo = [z1,z20]), z10)

Jz2_exact = convert(Matrix, ∂G_∂z2)
Jz2_num = FiniteDiff.finite_difference_jacobian(z2 -> G(μ0, x0; exo = [z10,z2]), z20)

print( maximum(abs, Jμ_num - Jμ_exact) )

print(maximum(abs, Jx_num - Jx_exact))

print(maximum(abs, Jz1_num - Jz1_exact) )

print(maximum(abs, Jz2_num - Jz2_exact))

# Plotting the differences between exact values calculated in ergordic.jl and numerical values calculated by FiniteDiff
using Plots

pl1 = spy(abs.(Jx_num).>1e-10, title="Numerical")
pl2 = spy(abs.(Jx_exact).>1e-10, title="Exact")
pl3 = spy(abs.(Jx_exact - Jx_num).>1e-10, title="Diff")
plot(pl1,pl2,pl3)

p1 = plot([Jz1_num,Jz1_exact], label=["Jz1_num" "Jz1_exact"])
p2 = plot([Jz1_num[1:50],Jz1_exact[1:50]], label=["Jz1_num_zoom" "Jz1_exact_zoom"])
p3 = plot((Jz1_num-Jz1_exact), label = "diff")
p4 = scatter((Jz1_num-Jz1_exact)[abs.(Jz1_num-Jz1_exact).>1e-8], label = "diff>1E-8", linestyle = :dot)
plot(p1,p2,p3,p4, layout = (4,1))


p1 = plot([Jz2_num,Jz2_exact], label=["Jz2_num" "Jz2_exact"])
p2 = plot([Jz2_num[31:60],Jz2_exact[31:60]], label=["num_zoom" "exact_zoom"])
p3 = plot((Jz2_num-Jz2_exact), label = "diff")
p4 = scatter((Jz2_num-Jz2_exact)[abs.(Jz2_num-Jz2_exact).>1e-8], label = "diff>1E-8", linestyle = :dot)
plot(p1,p2,p3,p4, layout = (4,1))


# Comparing the ergodic distributions when there is smoothing and when there is not

using PlotlyJS
using DataFrames

## rbc markov chain model

model = yaml_import("examples/models/rbc_mc.yaml")
sol = Dolo.improved_time_iteration(model)

nodes = sol.dr.grid.endo.nodes

k = [nodes[i][1] for i in 1:length(nodes)]

μ_no_smoothing = Dolo.ergodic_distribution(model, sol; smooth=false)
μ_smoothing = Dolo.ergodic_distribution(model, sol;  smooth=true)

n = Int(length(μ_smoothing)/2)
df1 = DataFrame(k=vcat(k,k), μ=vcat(2*μ_smoothing[1:n], 2*μ_no_smoothing[1:n]), smooth=vcat(["smoothing" for i in 1:n], ["no smoothing" for i in 1:n]))
df2 = DataFrame(k=vcat(k,k), μ=vcat(2*μ_smoothing[n+1:2*n], 2*μ_no_smoothing[n+1:2*n]), smooth=vcat(["smoothing" for i in 1:n], ["no smoothing" for i in 1:n]))


fig = PlotlyJS.make_subplots(
    rows=2, cols=1,
    column_widths=[1.],
    row_heights=[0.5, 0.5]
)

PlotlyJS.add_trace!(
    fig,
    PlotlyJS.scatter(df1[1:n,:],  x=:k, y=:μ, name  = "smoothing"),
    row=1, col = 1
    )

PlotlyJS.add_trace!(
    fig,
    PlotlyJS.scatter(df1[n+1:2*n,:],  x=:k, y=:μ, name = "no smoothing"),
    row=1, col = 1
    )

PlotlyJS.add_trace!(
    fig,
    PlotlyJS.scatter(df2[1:n,:],  x=:k, y=:μ, name = "smoothing"),
    row=2, col = 1
    )

PlotlyJS.add_trace!(
    fig,
    PlotlyJS.scatter(df2[n+1:2*n,:],  x=:k, y=:μ, name = "no smoothing"),
    row=2, col = 1
    )

relayout!(
    fig,
    margin=attr(r=10, t=25, b=40, l=60),
    annotations=[
        attr(
            text="k",
            showarrow=false,
            xref="paper",
            yref="paper",
            x=0.5,
            y=-0.05),
        attr(
            text="k",
            showarrow=false,
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.55),
        attr(
            text="Exogeneous shock 2",
            showarrow=false,
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.45),
        attr(
            text="Exogeneous shock 1",
            showarrow=false,
            xref="paper",
            yref="paper",
            x=0.5,
            y=1.05),
        attr(
            text="μ",
            showarrow=false,
            xref="paper",
            yref="paper",
            x=-0.05,
            y=0.8),
        attr(
            text="μ",
            showarrow=false,
            xref="paper",
            yref="paper",
            x=-0.05,
            y=0.2)
    ]
)

fig






## Aiyagari model
model = yaml_import("examples/models/consumption_savings_iid.yaml")
sol = Dolo.time_iteration(model)

nodes = sol.dr.grid.endo.nodes

w = [nodes[i][1] for i in 1:length(nodes)]

μ_no_smoothing = Dolo.ergodic_distribution(model, sol; smooth=false)
μ_smoothing = Dolo.ergodic_distribution(model, sol;  smooth=true)


n = Int(length(μ_smoothing))

df = DataFrame(w=vcat(w,w), μ=vcat(μ_smoothing[1:n], μ_no_smoothing[1:n]), smooth=vcat(["smoothing" for i in 1:n], ["no smoothing" for i in 1:n]))

PlotlyJS.plot(
    PlotlyJS.scatter(df, x=:w, y=:μ, group=:smooth),
    Layout(
        xaxis_title="w",
        yaxis_title="μ"
        )
        )



