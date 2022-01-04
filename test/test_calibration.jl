@testset "Testing model calibration" begin

    path = joinpath(Dolo.pkg_path, "examples", "models")
    model = Dolo.yaml_import(joinpath(path, "rbc_iid.yaml"))

    calib_0 = copy( model.calibration.flat )

    # check one value
    beta = model.calibration.flat[:beta]
    delta = model.calibration.flat[:delta]
    rk = model.calibration.flat[:rk]
    k0 = model.calibration[:states][2]
    @test beta==0.99
    @test rk==1/beta-1+delta
    @test abs(k0-9.35498)<1e-5


    # check calibrations are still correct after modifications
    Dolo.set_calibration!(model, :beta, 0.98)
    beta = model.calibration.flat[:beta]
    delta = model.calibration.flat[:delta]
    rk = model.calibration.flat[:rk]
    k0 = model.calibration[:states][2]
    @test beta==0.98
    @test rk==1/beta-1+delta
    @test abs(k0-6.37024)<1e-4

    # can we change many parameters at the same time?
    Dolo.set_calibration!(model; beta=0.97, i=:(delta*k-0.01))
    rk = model.calibration.flat[:rk]
    beta = model.calibration.flat[:beta]
    delta = model.calibration.flat[:delta]
    k0 = model.calibration[:states][2]
    i0 = model.calibration[:controls][2]
    @test beta==0.97
    @test rk==1/beta-1+delta
    @test i0==delta*k0-0.01

    # put everything back to original values
    set_calibration!(model; beta=0.99, i=:(delta*k))
    # check one value
    beta = model.calibration.flat[:beta]
    delta = model.calibration.flat[:delta]
    rk = model.calibration.flat[:rk]
    k0 = model.calibration[:states][2]
    i0 = model.calibration[:controls][2]
    @test beta==0.99
    @test rk==1/beta-1+delta
    @test abs(k0-9.35498)<1e-5
    @test i0==delta*k0




end
