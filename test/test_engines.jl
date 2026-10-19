const KernelAbstractions = Dolo.KernelAbstractions

@testset "Engines" begin

    engines = [Dolo.Threaded(), Dolo.KA(KernelAbstractions.CPU())]

    # largest difference between two vectors of SVector/SMatrix
    maxdiff(a, b) = maximum(maximum(abs.(u-v)) for (u,v) in zip(a,b))
    maxdiff(a::Dolo.LL, b::Dolo.LL) = maximum(
        maximum(abs.(u.F_x-v.F_x)) for (A,B) in zip(a.D, b.D) for (u,v) in zip(A,B)
    )

    for fn in ["rbc_mc", "rbc_ar1", "consumption_savings_iid"]

        @testset "$fn" begin

            model = include("$(pkgdir(Dolo))/examples/ymodels/$fn.jl")
            dm = Dolo.discretize(model)
            wk = Dolo.newton_workspace(dm)
            (; x0, r0, J, L, φ) = wk

            # reference values, with the serial engine
            r = Dolo.F!(deepcopy(r0), dm, x0, φ)
            J1 = Dolo.dF_1!(deepcopy(J), dm, x0, φ)
            L2 = Dolo.dF_2!(deepcopy(L), dm, x0, φ)
            Lv = Dolo.mul!(deepcopy(r0), L2, x0)
            Ls = Dolo.mul!(deepcopy(L2), -2.0)
            Ld = Dolo.ldiv!(J1, deepcopy(L2))

            for engine in engines
                @testset "$(nameof(typeof(engine)))" begin
                    @test maxdiff(Dolo.F!(deepcopy(r0), dm, x0, φ, engine).data, r.data) == 0
                    @test maxdiff(Dolo.dF_1!(deepcopy(J), dm, x0, φ, engine).data, J1.data) == 0
                    @test maxdiff(Dolo.dF_2!(deepcopy(L), dm, x0, φ, engine), L2) == 0
                    @test maxdiff(Dolo.mul!(deepcopy(r0), L2, x0, engine).data, Lv.data) < 1e-14
                    @test maxdiff(Dolo.mul!(deepcopy(L2), -2.0, engine), Ls) == 0
                    @test maxdiff(Dolo.ldiv!(J1, deepcopy(L2), engine), Ld) == 0
                end
            end

            # pointwise and allocating versions agree
            @test maxdiff(Dolo.F(dm, x0, φ).data, r.data) == 0
            @test maxdiff(Dolo.dF_2(dm, x0, φ), L2) == 0

        end
    end

    @testset "dF_1 with a linearized decision rule" begin
        # same jacobian as with dual numbers through the decision rule
        for (fn, disc) in [("rbc_mc", Dict()), ("integration_B", Dict(:endo => [4, 4, 4]))]
            model = include("$(pkgdir(Dolo))/examples/ymodels/$fn.jl")
            dm = Dolo.discretize(model, disc)
            (; x0, φ) = Dolo.time_iteration_workspace(dm; interp_mode=:cubic)
            x0.data .= [x .* (1 .+ 0.05*sin.((1:length(x)) .* n)) for (n,x) in enumerate(x0.data)]
            Dolo.fit!(φ, x0)
            d = maximum(
                maximum(abs.(Dolo.dF_1(dm, s, x, φ) - Dolo.ForwardDiff.jacobian(u -> Dolo.F(dm, s, u, φ), x)))
                for (s, x) in zip(Dolo.enum(dm.grid), x0.data)
            )
            @test d < 1e-12
        end
    end

    @testset "Blocks stored by component (BlockVector)" begin
        # same operations as with a vector of tuples of blocks
        model = include("$(pkgdir(Dolo))/examples/ymodels/rbc_mc.jl")
        dm = Dolo.discretize(model)
        (; x0, J, L, φ) = Dolo.newton_workspace(dm; interp_mode=:cubic)
        x0.data .= [x .* (1 + 0.1*sin(n)) for (n,x) in enumerate(x0.data)]
        Dolo.fit!(φ, x0)
        Dolo.dF_1!(J, dm, x0, φ)
        Dolo.dF_2!(L, dm, x0, φ)
        B = Dolo.LL(L.grid, Dolo.BlockVector(L.D), deepcopy(L.φ))
        @test collect(B.D) == L.D
        @test Dolo.dF_2!(deepcopy(B), dm, x0, φ).D == L.D
        @test maxdiff(Dolo.mul!(Dolo.duplicate(x0), B, x0).data, Dolo.mul!(Dolo.duplicate(x0), L, x0).data) == 0
        @test collect(Dolo.mul!(deepcopy(B), -2.0).D) == Dolo.mul!(deepcopy(L), -2.0).D
        @test maxdiff(Dolo.ldiv!(J, deepcopy(B)), Dolo.ldiv!(J, deepcopy(L))) < 1e-14
    end

    @testset "fit! updates the values of the decision rule" begin
        model = include("$(pkgdir(Dolo))/examples/ymodels/rbc_mc.jl")
        dm = Dolo.discretize(model)
        (; x0, φ) = Dolo.time_iteration_workspace(dm)
        ψ = deepcopy(φ)      # does not share its values with x0
        x0.data .= [x .* (1 + 0.1*sin(n)) for (n,x) in enumerate(x0.data)]
        Dolo.fit!(ψ, x0)
        @test ψ.values.data == x0.data
        @test all(ψ(s) == x for (s,x) in zip(Dolo.enum(dm.grid), x0.data))
    end

    @test Dolo.select_engine(:serial, nothing) isa Dolo.Serial
    @test Dolo.select_engine(:threads, nothing) isa Dolo.Threaded
    @test Dolo.select_engine(Dolo.Threaded(), nothing) isa Dolo.Threaded
    @test_throws ErrorException Dolo.select_engine(:foo, nothing)

end
