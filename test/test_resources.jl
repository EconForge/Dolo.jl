@testset "Resources" begin

    root = pkgdir(Dolo)
    report(d) = sprint(show, MIME"text/plain"(), d)

    d = Dolo.resources()
    @test d.system.threads == Threads.nthreads()
    @test occursin("System", report(d))

    # memory extrapolated from small grids = memory of the workspace on the grid itself
    model = include("$root/examples/ymodels/integration_B.jl")
    dm = Dolo.discretize(model, Dict(:endo => [7, 6, 5]))
    d = Dolo.resources(dm)
    parts, _ = Dolo.workspace_parts(dm)
    T, m = first(d.precisions)
    @test T == Float64
    @test [b for (_, b) in m.solvers] == [b for (_, b) in Dolo.solver_bytes(parts)]
    @test d.model.N == 9*7*6*5
    @test m.future_nodes == 9
    # Float32 copy, with all computations in Float32
    T32, m32 = d.precisions[2]
    @test T32 == Float32
    @test all(==(Float32), values(m32.results))
    @test occursin("Float32 copy of the model: all computations in Float32", report(d))

    # a model without Float32 copy (methods defined for typeof(model))
    model = include("$root/examples/ymodels/rbc_iid.jl")
    d = Dolo.resources(model)
    @test d.precisions[2][2] isa Exception
    @test any(a -> occursin("No Float32 copy", a[2]), d.advice)
    @test occursin("rbc_iid", report(d))

end
