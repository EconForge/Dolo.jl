import Dolo
import Dolang



@testset "features" begin

    path = Dolo.pkg_path

    examples_path = joinpath(path, "examples", "models")
    readdir(examples_path)

    model = Dolo.yaml_import(joinpath(examples_path, "rbc_dtcc_iid.yaml"))
    ff = features(model)
    @test ff[:one_dimensional] == true
    @test ff[:one_state] == true
    @test ff[:one_control] == true
    @test ff[:nonstochastic_transitions] == true
    @test ff[:bounds_are_constant] == false

    model = Dolo.yaml_import(joinpath(examples_path, "neoclassical.yaml"))
    ff = features(model)
    @test ff[:one_dimensional] == true
    @test ff[:one_state] == true
    @test ff[:one_control] == true
    @test ff[:nonstochastic_transitions] == true
    @test ff[:bounds_are_constant] == false

    model = Dolo.yaml_import(joinpath(examples_path, "rbc_dtcc_mc.yaml"))
    ff = features(model)
    @test f[:one_dimensional] == true
    @test f[:one_state] == true
    @test f[:one_control] == false
    @test f[:nonstochastic_transitions] == true
    @test f[:bounds_are_constant] == false

end
