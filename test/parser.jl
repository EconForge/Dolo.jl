@testset "Dolo._parse" begin
    @testset "numbers" begin
        for T in (Float16, Float32, Float64, Int8, Int16, Int32, Int64)
            x = rand(T)
            @test Dolo._parse(x) == x
        end
    end

    @testset "symbols" begin
        for i=1:10
            s = gensym()
            @test Dolo._parse(s) == symbol(string(s), "_")
        end
    end

    @testset "x_(shift_Integer)" begin
        for i=1:10, T in (Int8, Int16, Int32, Int64)
            @test Dolo._parse(string("x(", T(i), ")")) == symbol("x__$(i)_")
            @test Dolo._parse(string("x(", T(-i), ")")) == symbol("x_m$(i)_")
        end

        # only add underscore to naems when shift is 0
        @test Dolo._parse("x(0)") == :x_
    end

    @testset "other function calls" begin
        @testset "one argument" begin
            @test Dolo._parse("sin(x)") == :(sin(x_))
            @test Dolo._parse("sin(x(-1))") == :(sin(x_m1_))
            @test Dolo._parse("foobar(x(2))") == :(foobar(x__2_))
        end

        @testset "two arguments" begin
            @test Dolo._parse("dot(x, y(1))") == :(dot(x_, y__1_))
            @test Dolo._parse("plot(x(-1), y)") == :(plot(x_m1_, y_))
            @test Dolo._parse("bingbong(x(2), y)") == :(bingbong(x__2_, y_))
        end

        @testset "more args" begin
            for i=3:10
                ex = Expr(:call, :(+), [:(x($j)) for j in 1:i]...)
                want = Expr(:call, :(+), [symbol("x__", j, "_") for j in 1:i]...)
                @test Dolo._parse(ex) == want
            end
        end

        @testset "throws errors when unsupported" begin
            @test_throws ErrorException Dolo._parse("x+y | i <= 100")
        end
    end

    @testset "Expr(:(=), ...)" begin
        @testset "without targets" begin
            @test Dolo._parse(:(x = y)) == :(y_ - x_)
        end

        @testset "with targets" begin
            @test Dolo._parse(:(x = log(y(-1))); targets=[:x]) == :(x_ = log(y_m1_))
            @test_throws ErrorException Dolo._parse(:(x = y); targets=[:y])
        end
    end
end

type MockSymbolic <: Dolo.ASM
    symbols
    equations
end

MockSymbolic(d) = MockSymbolic(d, nothing)
model_spec(::MockSymbolic) = :dtcscc
@testset "compiler" begin

    @testset "_param_block" begin
        sm = MockSymbolic(Dict(:parameters => [:a, :b, :foobar]))
        have = Dolo._param_block(sm)
        @test have.head == :block
        @test have.args[1] == :(@inbounds a_ = p[1])
        @test have.args[2] == :(@inbounds b_ = p[2])
        @test have.args[3] == :(@inbounds foobar_ = p[3])
    end

    @testset "_single_arg_block" begin
        sm = MockSymbolic(Dict(:states => [:z, :k]))

        @testset "no shift" begin
            have = Dolo._single_arg_block(sm, :s, :states, 0)
            @test have.head == :block
            @test have.args[1] == :(@inbounds z_ = s[1])
            @test have.args[2] == :(@inbounds k_ = s[2])
        end

        @testset "positive shift" begin
            have = Dolo._single_arg_block(sm, :s, :states, 1)
            @test have.head == :block
            @test have.args[1] == :(@inbounds z__1_ = s[1])
            @test have.args[2] == :(@inbounds k__1_ = s[2])
        end

        @testset "negative shift" begin
            have = Dolo._single_arg_block(sm, :s, :states, -1)
            @test have.head == :block
            @test have.args[1] == :(@inbounds z_m1_ = s[1])
            @test have.args[2] == :(@inbounds k_m1_ = s[2])
        end
    end

    @testset "_assign_single_el" begin
        @test Dolo._assign_single_el(:out, Inf, 1) == :(out[1] = $Inf)
        @test Dolo._assign_single_el(:z__m_1, 0, :foo) == :(z__m_1[foo]=0)
    end

    @testset "_main_body_block" begin
        sm = MockSymbolic(Dict())
        @testset "with targets" begin
            targets = [:x, :y, :z]
            exprs = [:(x = sin(42.0)), :(y = x(-1) + cos(42.0)), :(z = tan(x))]
            have = Dolo._main_body_block(sm, targets, exprs)
            @test have.head == :block
            @test have.args[1] == :(out[1] = x_ = sin(42.0))
            @test have.args[2] == :(out[2] = y_ = x_m1_ + cos(42.0))
            @test have.args[3] == :(out[3] = z_ = tan(x_))
        end

        @testset "throws proper errors" begin
            # can't have shift on lhs
            tar = [:z]
            ex = [:(z(1) = 10.0)]
            @test_throws ErrorException Dolo._main_body_block(sm, tar, ex)

            # can't have two symbols on lhs
            ex = [:(z + y = 10.0)]
            @test_throws ErrorException Dolo._main_body_block(sm, tar, ex)

            # lhs not in targets
            tar = [:y]
            ex = [:(z = 10.0)]
            @test_throws ErrorException Dolo._main_body_block(sm, tar, ex)
        end

        @testset "without targets" begin
            targets = Symbol[]
            exprs = [:(x = sin(42.0)), :(y = x(-1)+cos(42.0)), :(z+y = tan(x))]
            have = Dolo._main_body_block(sm, targets, exprs)
            @test have.head == :block
            @test have.args[1] == :(out[1] = sin(42.0) - x_)
            @test have.args[2] == :(out[2] = x_m1_ + cos(42.0) - y_)
            @test have.args[3] == :(out[3] = tan(x_) - (z_ + y_))
        end
    end

    # the next two testsets use this data
    exprs = [parse("y = z*k^alpha*n^(1-alpha)"),
             parse("c = y - i"),
             parse("rk = alpha*y/k"),
             parse("w = (1-alpha)*y/n")]
    sm = MockSymbolic(Dict(:auxiliaries => [:y, :rk, :w, :c],
                           :states => [:z, :k],
                           :controls => [:i, :n],
                           :parameters => [:fizz, :buzz, :alpha]),
                      Dict(:auxiliary => exprs, :value => Expr[]))

    @testset "_aux_block" begin
        @testset "without shift" begin
            have = Dolo._aux_block(sm, 0)
            @test have.head == :block
            @test have.args[1] == :(y_ = z_*k_^alpha_*n_^(1-alpha_))
            @test have.args[2] == :(c_ = y_ - i_)
            @test have.args[3] == :(rk_ = alpha_*y_/k_)
            @test have.args[4] == :(w_ = (1-alpha_)*y_/n_)
        end

        @testset "positive shift" begin
            have = Dolo._aux_block(sm, 1)
            @test have.head == :block
            @test have.args[1] == :(y__1_ = z__1_*k__1_^alpha_*n__1_^(1-alpha_))
            @test have.args[2] == :(c__1_ = y__1_ - i__1_)

            # notice below that we have k and rk, but that the single shift was
            # properly applied one time to each of them. Had we not used regex
            # in _aux_block and just searched for the state, then control, then
            # auxiliary name we would have gotten rk(1)(1) instead of rk(1)
            @test have.args[3] == :(rk__1_ = alpha_*y__1_/k__1_)
            @test have.args[4] == :(w__1_ = (1-alpha_)*y__1_/n__1_)
        end

        @testset "negative shift" begin
            have = Dolo._aux_block(sm, -2)
            @test have.head == :block
            @test have.args[1] == :(y_m2_ = z_m2_*k_m2_^alpha_*n_m2_^(1-alpha_))
            @test have.args[2] == :(c_m2_ = y_m2_ - i_m2_)
            @test have.args[3] == :(rk_m2_ = alpha_*y_m2_/k_m2_)
            @test have.args[4] == :(w_m2_ = (1-alpha_)*y_m2_/n_m2_)
        end
    end

    @testset "compile_equation" begin
        # can't compile this type of equation, so calling throws an error.
        # it is still a dolo functor, however
        obj1 = let
            eval(compile_equation(sm, :value))
        end
        @test_throws ErrorException obj1()
        @test super(typeof(obj1)) == Dolo.AbstractDoloFunctor

        obj2 = let
            eval(compile_equation(sm, :auxiliary))
        end
        @test super(typeof(obj2)) == Dolo.AbstractDoloFunctor

        # now see if the functions were compiled properly
        z = 1.0
        k = 0.25
        i = 100.0
        n = 10.0
        alpha = 0.33
        s = [z, k]
        x = [i, n]
        p = [0.0, 0.0, alpha]

        y = z*k^alpha*n^(1-alpha)
        c = y-i
        rk = alpha*y/k
        w = (1-alpha)*y/n

        want = [y, c, rk, w]
        out = zeros(4)

        @test_throws MethodError obj2()
        @test_throws MethodError obj2(s)
        @test_throws MethodError obj2(s, x)

        # test allocating version
        @test @inferred(obj2(s, x, p)) == want

        # test non-allocating version
        @inferred obj2(s, x, p, out)
        @test out == want
    end

end
