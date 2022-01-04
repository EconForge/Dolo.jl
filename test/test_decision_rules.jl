@testset "DR types" begin
    grid_smol = Dolo.SmolyakGrid{2}([0.0, 0.0], [1.0, 1.0,], [3,3])
    grid_exo = Dolo.EmptyGrid{2}()
    grid_exo2 = Dolo.UnstructuredGrid{2}(rand(1, 2))
    grid_rand = Dolo.RandomGrid{2}([0.0, 0.0], [1.0, 1.0,], 50)
    pts = rand(200, 2)
    pts_vec = copy(reshape(reinterpret(Dolo.Point{2}, vec(copy(pts'))), (size(pts, 1),)))
    test_sets = [
        (grid_smol, Dolo.SmolyakDR, Dolo.Smolyak),
        (grid_rand, Dolo.CompletePolyDR, Dolo.CompletePolnomial{3}),
        (grid_rand, Dolo.CompletePolyDR, Dolo.CompletePolnomial),
        (grid_rand, Dolo.CompletePolyDR, Dolo.CompletePolnomial{2}),
        (grid_smol, Dolo.CompletePolyDR, Dolo.CompletePolnomial{3}),
        (grid_smol, Dolo.CompletePolyDR, Dolo.CompletePolnomial),
        (grid_smol, Dolo.CompletePolyDR, Dolo.CompletePolnomial{2}),
        (grid_exo2, Dolo.CompletePolyDR, Dolo.CompletePolnomial{3}),
        (grid_exo2, Dolo.CompletePolyDR, Dolo.CompletePolnomial),
        (grid_exo2, Dolo.CompletePolyDR, Dolo.CompletePolnomial{2}),
    ]

    @testset for (grid_endo, dr_name, drT) in test_sets

       g = Dolo.nodes(Matrix, grid_endo)
       sg_vals = [sin.(g[:, 1]) cos.(g[:, 2])]
       sg_vals_vec = copy(reshape(reinterpret(Dolo.Point{2}, vec(copy(sg_vals'))), (size(g,1),)))

        dr1 = dr_name(grid_exo, grid_endo, Val{2})
        dr2 = Dolo.DecisionRule(grid_exo, grid_endo, Val{2}, drT)
        for dr in [dr1, dr2]
            Dolo.set_values!(dr, [sg_vals])
            out = Dolo.evaluate(dr, pts)
            @test maximum(abs, out - Dolo.evaluate(dr, pts_vec)) < 1e-14  ### ?
            Dolo.set_values!(dr, [sg_vals_vec])
            @test maximum(abs, out - Dolo.evaluate(dr, pts_vec)) < 1e-14
            @test maximum(abs, out - Dolo.evaluate(dr, pts)) < 1e-14
            @test_throws ErrorException Dolo.set_values!(dr, [sg_vals_vec, sg_vals_vec])
        end

        dr1 = dr_name(grid_exo2, grid_endo, Val{2})
        dr2 = dr_name(grid_exo2, grid_endo, Val{2})
        for dr in [dr1, dr2]
            Dolo.set_values!(dr, [sg_vals])
            out = Dolo.evaluate(dr, 1, pts)
            @test maximum(abs, out - Dolo.evaluate(dr, 1, pts_vec)) < 1e-14
            Dolo.set_values!(dr, [sg_vals_vec])
            @test maximum(abs, out - Dolo.evaluate(dr, 1, pts_vec)) < 1e-14
            @test maximum(abs, out - Dolo.evaluate(dr, 1, pts)) < 1e-14
            @test_throws MethodError Dolo.evaluate(dr, pts)
            @test_throws MethodError Dolo.evaluate(dr, pts_vec)
            @test_throws MethodError Dolo.evaluate(dr, 1, pts_vec[1])
            @test_throws BoundsError Dolo.evaluate(dr, 2, pts_vec)
            @test_throws ErrorException Dolo.set_values!(dr, [sg_vals_vec, sg_vals_vec])
        end
    end
end
