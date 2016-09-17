
function import_model(fname)
    yaml_import(fname)
    return true
end

# @test import_model("/home/pablo/Programming/econforge/dolo/examples/models/rbc.yaml")

@testset "Testing yaml_import" begin
    # @test import_model("https://raw.githubusercontent.com/EconForge/dolo/144965224f432c9f467f0e667bc0cc4d77caf629/examples/models/rbc.yaml")
    @test import_model(Pkg.dir("Dolo", "examples", "models", "rbc_dtcc.yaml"))
    # TODO: re-enable once we can build objects of type MarkovChain
    # @test import_model("https://raw.githubusercontent.com/EconForge/dolo/144965224f432c9f467f0e667bc0cc4d77caf629/examples/models/sudden_stop.yaml")
end
