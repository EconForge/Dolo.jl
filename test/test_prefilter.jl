@testset "Spline prefilter" begin

    splines = Dolo.splines
    SVector = Dolo.SVector

    # random values with padding, as in a SplineInterpolator
    function padded(T, n...)
        θ = zeros(T, (k+2 for k in n)...)
        θ[(2:k+1 for k in n)...] .= rand(T, n...)
        return θ
    end

    # chunked version (here with KernelAbstractions on CPU) == sequential version
    @testset "chunked ≈ sequential: $T, size $n" for T in (Float64, Float32, SVector{2,Float64}),
            n in ((5,), (100,), (1000,), (300, 7), (20, 150), (6, 140, 5))
        θ = padded(T, n...)
        seq = splines.prefilter_sequential!(copy(θ))
        chk = splines.prefilter_chunked!(copy(θ), similar(θ))
        tol = 50*eps(eltype(eltype(θ)))
        @test maximum(maximum(abs.(u-v)) for (u,v) in zip(seq, chk)) < tol
    end

    # in 1D, the sequential version is prefilter_2!
    θ = padded(Float64, 1000)
    @test splines.prefilter_sequential!(copy(θ)) == (c = copy(θ); splines.prefilter_2!(c); c)

end
