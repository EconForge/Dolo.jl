
function import_model(fname)
    yaml_import(fname)
    return true
end

# @test import_model("/home/pablo/Programming/econforge/dolo/examples/models/rbc.yaml")

@testset "Testing yaml_import" begin
    @test import_model("https://raw.githubusercontent.com/EconForge/dolo/c8bd2e3f2f5402f687beb7949d49deefda6a5fc6/examples/models/rbc.yaml")
    @test import_model("https://raw.githubusercontent.com/EconForge/dolo/c8bd2e3f2f5402f687beb7949d49deefda6a5fc6/examples/models/sudden_stop.yaml")
end
