
function import_model(fname)
    yaml_import(fname)
    return true
end

# @test import_model("/home/pablo/Programming/econforge/dolo/examples/models/rbc.yaml")

@testset "Testing yaml_import" begin
    # @test import_model("https://raw.githubusercontent.com/EconForge/dolo/144965224f432c9f467f0e667bc0cc4d77caf629/examples/models/rbc.yaml")
    path = Pkg.dir("Dolo", "examples", "models")
    files = [f for f in readdir(path) if contains(f, "yaml")]
    for fname in files
        print("Importing ", fname,"\n")
        model = yaml_import(joinpath(path, fname))
        res = residuals(model)
        if !contains(fname,"sudden_stop")
            @test (maximum(res[:transition]) + maximum(res[:arbitrage]) )<1e-5
        else
            @test (maximum(res[:transition]) + maximum(res[:arbitrage])-0.215)<1e-5
        end
    end

    # TODO: re-enable once we can build objects of type MarkovChain
    # @test import_model("https://raw.githubusercontent.com/EconForge/dolo/144965224f432c9f467f0e667bc0cc4d77caf629/examples/models/sudden_stop.yaml")
end
