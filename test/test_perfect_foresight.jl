import Dolo

path = Dolo.pkg_path

@testset "testing perfect_foresight" begin

    joinpath(path, "examples", "models", "rbc_dtcc_ar1.yaml");

    # Test with two shocks
    fn = joinpath(path, "examples", "models", "rbc_catastrophe.yaml");

    model = Dolo.yaml_import(fn);

    # Provide the shocks
    # If not provided will be equal to zero by default
    # Vector sizes for each shock can be different
    # The order the shocks are defined are irrelevant: pf fuction finds the index for each shock

    exo = Dict(:z =>  [0.00, 0.01, 0.02, 0.03, 0.04], :xi =>  [0, 0, 0, 0, 0, 0])
    # or
    exo = Dict(:xi =>  [0, 0, 0, 0, 0, 0] , :z =>  [0.00, 0.01, 0.02, 0.03, 0.04])
    # or
    exo = Dict(:z =>  [0.00, 0.01, 0.02, 0.03, 0.04])

    @time df = Dolo.perfect_foresight(model, exo, T=20)

end
