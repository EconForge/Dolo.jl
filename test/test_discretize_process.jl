import Dolo

path = Dolo.pkg_path

@testset "testing discretizing functions" begin

    fn = joinpath(path, "LAMP_2s.yaml")
    model = Dolo.yaml_import(fn)

    Dolo.discretize(model.exogenous)
    exo = Dolo.discretize(Dolo.DiscretizedProcess, model.exogenous)
    # @time sol = Dolo.time_iteration(model, exo; maxit=20, verbose=false, maxit=10000)


    @test true

end
