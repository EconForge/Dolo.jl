import Dolo

model = Dolo.yaml_import("examples/models/rbc_dtcc_mc.yaml")

# @time sol = Dolo.time_iteration(model, maxit=1000)


@time ret = Dolo.improved_time_iteration(model, complementarities=true)


# @time ddr = Dolo.improved_time_iteration(model, complementarities=true)

vvvv = Dolo.CachedDecisionRule{Dolo.CubicDR{Dolo.EmptyGrid,Dolo.CartesianGrid{2},2,2},Dolo.DiscretizedIIDProcess}




(::Int64, ::Int64, ::Array{SVector{2,Float64},1})
