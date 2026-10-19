const mul! = Dolo.mul!
const norm = Dolo.norm
eye(n) = [i==j ? 1.0 : 0.0 for i=1:n, j=1:n]

@testset "Linear solvers" begin

    @testset "GMRES (dense)" begin
        n = 40
        A = eye(n) + 0.3*sin.((1:n) .* (1:n)')/sqrt(n)
        b = cos.(1:n)
        A! = (w,v) -> mul!(w, A, v)
        wk = Dolo.gmres_workspace(Float64, n, n)

        x = zeros(n)
        k, conv = Dolo.gmres!(x, A!, b, wk; tol=1e-12)
        @test conv
        @test k <= n
        @test x ≈ A\b

        # Krylov space too small: not converged, but residual decreased
        wk5 = Dolo.gmres_workspace(Float64, n, 5)
        x = zeros(n)
        k, conv = Dolo.gmres!(x, A!, b, wk5; tol=1e-12)
        @test k == 5
        @test norm(A*x-b) < norm(b)

        # lucky breakdown: exact solution after one iteration
        x = zeros(n)
        k, conv = Dolo.gmres!(x, (w,v)->(w .= 2 .* v), b, wk; tol=1e-12)
        @test k == 1
        @test conv
        @test x ≈ b/2

        # initial guess is used
        x = A\b
        k, conv = Dolo.gmres!(x, A!, b, wk; tol=1e-10)
        @test k == 0
    end

    @testset "(I-L) dx = r" begin

        model = include("$(pkgdir(Dolo))/examples/ymodels/rbc_mc.jl")
        dm = Dolo.discretize(model)
        wk = Dolo.time_iteration_workspace(dm; accelerate=true, linsolve=:gmres, krylov_dim=100)
        (;x0, r0, J, L, φ) = wk
        Dolo.dF_1!(J, dm, x0, φ)
        Dolo.dF_2!(L, dm, x0, φ)
        Dolo.mul!(L, -1.0)
        Dolo.ldiv!(J, L)

        M = convert(Matrix, L)
        r0.data .= x0.data
        r = copy(Dolo.ravel(r0))
        sol = (eye(size(M,1)) - M) \ r

        for linsolve in (:neumann, :gmres)
            dx = Dolo.duplicate(x0)
            k = Dolo.solve_I_minus_L!(dx, L, r0, wk, Dolo.Serial(); linsolve=linsolve, tol=1e-12)
            @test k > 0
            @test Dolo.ravel(dx) ≈ sol rtol=1e-8
        end

    end

end
