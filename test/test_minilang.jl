@testset "minilang.jl tests" begin
    mc = Dolo.ModelCalibration(
        OrderedDict{Symbol,Real}(:sigz => 0.1, :sigy => 0.2),
        OrderedDict(:exog => [:sigz, :sigy])
    )
    good_data = Dict(
        :tag => :Normal,
        :Sigma => [[:(sigz^2), 0.0], [0.0, :(sigy^2)]]
    )

    bad_data = Dict(
        :tag => :Normal,
        :Sigma => [[:(sigz^2)], [:(sigy^2)]]
    )

    @test_throws DimensionMismatch Dolo._build_exogenous_entry(bad_data, mc)

    want = Dolo.MvNormal([0.0, 0.0], [0.01 0.0; 0.0 0.04])
    @test begin
        have = Dolo._build_exogenous_entry(good_data, mc)
        have.mu ≈ want.mu && want.Sigma ≈ have.Sigma
    end
end
