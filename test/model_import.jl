@testset "Testing yaml_import" begin
    path = joinpath(Dolo.pkg_path, "examples", "models")
    files = [f for f in readdir(path) if occursin("yaml", f)]
    for fname in files
        println("Importing ", fname)
        model = yaml_import(joinpath(path, fname))
        res = residuals(model)
        if !occursin("sudden_stop", fname)
            @test (maximum(res[:transition]) + maximum(res[:arbitrage]) )<1e-5
        else
            @test (maximum(res[:transition]) + maximum(res[:arbitrage])-0.215)<1e-5
        end
    end
end
