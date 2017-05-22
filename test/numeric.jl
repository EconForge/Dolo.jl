@testset "TaylorExpansion" begin
    s0 = [0.2, 0.4, 1.1]
    x0 = [1.2, 0.9]

    N = 1000
    points = rand(N, 3)

    X_s = rand(2, 3)
    X_ss = rand(2, 3, 3)
    X_sss = rand(2, 3, 3, 3)
    dr1 = TaylorExpansion(s0, x0, X_s)
    dr2 = TaylorExpansion(s0, x0, X_s, X_ss)
    dr3 = TaylorExpansion(s0, x0, X_s, X_ss, X_sss)

    out1 = dr1(points)
    out2 = dr2(points)
    out3 = dr3(points)

    out1_1d = dr1(vec(points[1, :]))
    out2_1d = dr2(vec(points[1, :]))
    out3_1d = dr3(vec(points[1, :]))

    @test maximum(abs, out1_1d - vec(out1[1, :])) == 0.0
    @test maximum(abs, out2_1d - vec(out2[1, :])) == 0.0
    @test maximum(abs, out3_1d - vec(out3[1, :])) == 0.0

    ds = points .- s0'
    verif1 = x0 .+ X_s*ds'
    @test maximum(abs, out1 - verif1') < 1e-12

    # TODO: write verification tests for 2d and 3d
end
