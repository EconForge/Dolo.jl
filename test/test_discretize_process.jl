import Dolo

path = Dolo.pkg_path

@testset "testing discretizing functions" begin

    rho_z = 0.9
    rho_b = 0.9
    sig_z = 0.001
    sig_b = 0.001

    proc = Dolo.ProductProcess(
        Dolo.VAR1(rho_z, ones(1,1)*sig_z^2),
        Dolo.VAR1(rho_b, ones(1,1)*sig_b^2)
    )

    Dolo.discretize(proc)
    exo = Dolo.discretize(Dolo.GDP, proc)

    @test true

end
