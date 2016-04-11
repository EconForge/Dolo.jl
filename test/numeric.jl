@testset "TaylorExpansion" begin
    s0 = [0.2, 0.4, 1.1]
    x0 = [1.2, 0.9]

    N = 1000
    points = rand(3, N)

    X_s = rand(2, 3)
    X_ss = rand(2, 3, 3)
    X_sss = rand(2, 3, 3, 3)
    dr1 = TaylorExpansion(s0, x0, X_s)
    dr2 = TaylorExpansion(s0, x0, X_s, X_ss)
    dr3 = TaylorExpansion(s0, x0, X_s, X_ss, X_sss)

    out1 = dr1(points)
    out2 = dr2(points)
    out3 = dr3(points)

    out1_1d = dr1(points[:, 1])
    out2_1d = dr2(points[:, 1])
    out3_1d = dr3(points[:, 1])

    @test maxabs(out1_1d - out1[:, 1]) == 0.0
    @test maxabs(out2_1d - out2[:, 1]) == 0.0
    @test maxabs(out3_1d - out3[:, 1]) == 0.0

    ds = points .- s0
    verif1 = x0 .+ X_s*ds
    @test maxabs(out1 - verif1) < 1e-12

    # TODO: write verification tests for 2d and 3d
end
