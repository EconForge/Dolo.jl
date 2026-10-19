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
