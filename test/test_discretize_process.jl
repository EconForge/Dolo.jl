import Dolo

path = Dolo.pkg_path

@testset "testing discretizing functions" begin

    rho_z = 0.9
    rho_b = 0.9
    sig_z = 0.001
    sig_b = 0.001

    proc = Dolo.ProductProcess(
        Dolo.VAR1(rho_z, eye(1)*sig_z^2),
        Dolo.VAR1(rho_b, eye(1)*sig_b^2))

    Dolo.discretize(proc)
    exo = Dolo.discretize(Dolo.DiscretizedProcess, proc)

    @test true

end
