@testset "Testing yaml_import" begin
    path = joinpath(Dolo.pkg_path, "examples", "models")
    files = [f for f in readdir(path) if contains(f, "yaml")]
    for fname in files
        println("Importing ", fname)
        model = yaml_import(joinpath(path, fname))
        res = residuals(model)
        if !contains(fname,"sudden_stop")
            @test (maximum(res[:transition]) + maximum(res[:arbitrage]) )<1e-5
        else
            @test (maximum(res[:transition]) + maximum(res[:arbitrage])-0.215)<1e-5
        end
    end
end
