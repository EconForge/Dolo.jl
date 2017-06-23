import Dolo

@testset "features" begin

    path = Dolo.pkg_path

    examples_path = joinpath(path, "examples", "models")
    readdir(examples_path)

    model = Dolo.yaml_import(joinpath(examples_path, "rbc_dtcc_iid.yaml"))
    ff = Dolo.features(model)
    @test ff[:one_dimensional] == false
    @test ff[:one_state] == false
    @test ff[:one_control] == false
    @test ff[:nonstochastic_transitions] == false
    @test ff[:bounds_are_constant] == true

    model = Dolo.yaml_import(joinpath(examples_path, "neoclassical.yaml"))
    ff = Dolo.features(model)
    @test ff[:one_dimensional] == true
    @test ff[:one_state] == true
    @test ff[:one_control] == true
    @test ff[:nonstochastic_transitions] == true
    @test ff[:bounds_are_constant] == true

    model = Dolo.yaml_import(joinpath(examples_path, "rbc_dtcc_mc.yaml"))
    ff = Dolo.features(model)
    @test ff[:one_dimensional] == true
    @test ff[:one_state] == true
    @test ff[:one_control] == false
    @test ff[:nonstochastic_transitions] == true
    @test ff[:bounds_are_constant] == true

end
