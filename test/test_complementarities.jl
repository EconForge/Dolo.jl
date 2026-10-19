@testset "Complementarities" begin

    ⫫ = Dolo.:⫫
    frule = Dolo.ChainRulesCore.frule

    @test ⫫(0.0, 0.0) == 0.0
    # (u+v-√(u²+v²))/2 → u/2 when v → ∞
    @test ⫫(1.0, Inf) ≈ 0.5 atol=1e-4
    @test ⫫(Inf, 1.0) ≈ 0.5 atol=1e-4

    # frule agrees with finite differences away from the kink
    for (u,v) in ((0.3, -0.2), (1.0, 2.0), (-0.5, 0.1))
        h = 1e-6
        y, dy_u = frule((nothing, 1.0, 0.0), ⫫, u, v)
        _, dy_v = frule((nothing, 0.0, 1.0), ⫫, u, v)
        @test y == ⫫(u,v)
        @test dy_u ≈ (⫫(u+h,v)-⫫(u-h,v))/(2h) atol=1e-6
        @test dy_v ≈ (⫫(u,v+h)-⫫(u,v-h))/(2h) atol=1e-6
    end

    # well defined at the origin and for infinite bounds
    @test frule((nothing, 1.0, 1.0), ⫫, 0.0, 0.0)[2] == 0.0
    @test frule((nothing, 1.0, 1.0), ⫫, 1.0, Inf)[2] ≈ 0.5 atol=1e-4

end
