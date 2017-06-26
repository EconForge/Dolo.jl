import Dolo

path = Dolo.pkg_path

@testset "testing perfect_foresight" begin

    fn = joinpath(path, "examples", "models", "rbc_dtcc_ar1.yaml")

    model = Dolo.yaml_import(fn)

    # we must define the series for exogenous shocks
    n_e = length(model.symbols[:exogenous])

    T_e = 5
    exo = zeros(T_e, n_e)
    exo[1,:] = [0.00]  # this is used to determine initial steady-state of the model
    exo[2,:] = [0.01]
    exo[3,:] = [0.02]
    exo[4,:] = [0.03]
    exo[5,:] = [0.04]  # this is used to determine final steady-state

    @time df = Dolo.perfect_foresight(model, exo, T=20)

    @test true

end
