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
                ex = Expr(:call, :my_func, [:(x($j)) for j in 1:i]...)
                want = Expr(:call, :my_func, [symbol("x__", j, "_") for j in 1:i]...)
                @test Dolo._parse(ex) == want
            end
        end

        @testset "arithmetic" begin
            @test Dolo._parse(:(a(1) + b + c(2) + d(-1))) == :(((a__1_ .+ b_) .+ c__2_) .+ d_m1_)
            @test Dolo._parse(:(a(1) * b * c(2) * d(-1))) == :(((a__1_ .* b_) .* c__2_) .* d_m1_)
            @test Dolo._parse(:(a(1) - b - c(2) - d(-1))) == :(((a__1_ .- b_) .- c__2_) .- d_m1_)
            @test Dolo._parse(:(a(1) ^ b)) == :(a__1_ .^ b_)
        end

        @testset "throws errors when unsupported" begin
            @test_throws ErrorException Dolo._parse("x+y || i <= 100")
        end
    end

    @testset "Expr(:(=), ...)" begin
        @testset "without targets" begin
            @test Dolo._parse(:(x = y)) == :(y_ .- x_)
        end

        @testset "with targets" begin
            @test Dolo._parse(:(x = log(y(-1))); targets=[:x]) == :(x_ = log(y_m1_))
            @test_throws ErrorException Dolo._parse(:(x = y); targets=[:y])
        end
    end
end

type MockSymbolic{ID,kind} <: Dolo.ASM{ID,kind}
    symbols
    equations
end

MockSymbolic(d, eq) = MockSymbolic{:foobar,:dtcscc}(d, eq)
MockSymbolic(d) = MockSymbolic(d, nothing)
function Dolo.DTCSCCModel{id}(sm::MockSymbolic{id})
    calib = ModelCalibration(FlatCalibration(Dolo.OrderedDict()),
                             GroupedCalibration(Dict()),
                             Dict(),
                             OrderedDict())
    d = Dict{Symbol,Any}()
    DTCSCCModel{id}(sm, calib, d, d, :dtcscc, "foobar", "boo.yaml")
end

@testset "compiler" begin

    @testset "_param_block" begin
        sm = MockSymbolic(Dict(:parameters => [:a, :b, :foobar]))
        have = Dolo._param_block(sm)
        @test have.head == :block
        @test have.args[1] == :(a_ = _unpack_var(p, 1))
        @test have.args[2] == :(b_ = _unpack_var(p, 2))
        @test have.args[3] == :(foobar_ = _unpack_var(p, 3))
    end

    @testset "_single_arg_block" begin
        sm = MockSymbolic(Dict(:states => [:z, :k]))

        @testset "no shift" begin
            have = Dolo._single_arg_block(sm, :s, :states, 0)
            @test have.head == :block
            @test have.args[1] == :(z_ = _unpack_var(s, 1))
            @test have.args[2] == :(k_ = _unpack_var(s, 2))
        end

        @testset "positive shift" begin
            have = Dolo._single_arg_block(sm, :s, :states, 1)
            @test have.head == :block
            @test have.args[1] == :(z__1_ = _unpack_var(s, 1))
            @test have.args[2] == :(k__1_ = _unpack_var(s, 2))
        end

        @testset "negative shift" begin
            have = Dolo._single_arg_block(sm, :s, :states, -1)
            @test have.head == :block
            @test have.args[1] == :(z_m1_ = _unpack_var(s, 1))
            @test have.args[2] == :(k_m1_ = _unpack_var(s, 2))
        end
    end

    @testset "_assign_single_el" begin
        @test Dolo._assign_var_expr(:out, Inf, 1) == :(_assign_var(out, $Inf, 1))
        @test Dolo._assign_var_expr(:z__m_1, 0, :foo) == :(_assign_var(z__m_1, 0, foo))
    end

    @testset "_main_body_block" begin
        sm = MockSymbolic(Dict())
        @testset "with targets" begin
            targets = [:x, :y, :z]
            exprs = [:(x = sin(42.0)), :(y = x(-1) + cos(42.0)), :(z = tan(x))]
            have = Dolo._main_body_block(sm, targets, exprs)
            @test have.head == :block
            @test have.args[1] == :(x_ = sin(42.0))
            @test have.args[2] == :(y_ = x_m1_ .+ cos(42.0))
            @test have.args[3] == :(z_ = tan(x_))
            @test have.args[4] == :(_assign_var(out, x_, 1))
            @test have.args[5] == :(_assign_var(out, y_, 2))
            @test have.args[6] == :(_assign_var(out, z_, 3))
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
            @test have.args[1] == :(_assign_var(out, sin(42.0) .- x_, 1))
            @test have.args[2] == :(_assign_var(out, x_m1_ .+ cos(42.0) .- y_, 2))
            @test have.args[3] == :(_assign_var(out, tan(x_) .- (z_ .+ y_), 3))
        end
    end

    # the next two testsets use this data
    # NOTE: I know that transition specifies equations in the wrong order.
    # NOTE: I know that the functions for expecation have too many variables
    #       these are part of the tests below
    exprs = [parse("y = z*k^alpha*n^(1-alpha)"),
             parse("c = y - i"),
             parse("rk = alpha*y/k"),
             parse("w = (1-alpha)*y/n")]
    sm = MockSymbolic(Dict(:auxiliaries => [:y, :c, :rk, :w],
                           :states => [:z, :k],
                           :shocks => [:ɛz],
                           :controls => [:i, :n],
                           :expectations => [:Ez],
                           :parameters => [:fizz, :buzz, :alpha]),
                      Dict(:auxiliary => exprs, :value => Expr[],
                           :transition => [:(k = k(-1) + 1), :(z = z(-1))],
                           :expectation => [:(Ez = z(1)), :(Ek = k(1))]))

    m = DTCSCCModel(sm)

    @testset "_aux_block" begin
        @testset "without shift" begin
            have = Dolo._aux_block(sm, 0)
            @test have.head == :block
            @test have.args[1] == :(y_ = z_.*k_.^alpha_.*n_.^(1.-alpha_))
            @test have.args[2] == :(c_ = y_ .- i_)
            @test have.args[3] == :(rk_ = alpha_.*y_./k_)
            @test have.args[4] == :(w_ = (1.-alpha_).*y_./n_)
        end

        @testset "positive shift" begin
            have = Dolo._aux_block(sm, 1)
            @test have.head == :block
            @test have.args[1] == :(y__1_ = z__1_.*k__1_.^alpha_.*n__1_.^(1.-alpha_))
            @test have.args[2] == :(c__1_ = y__1_ .- i__1_)

            # notice below that we have k and rk, but that the single shift was
            # properly applied one time to each of them. Had we not used regex
            # in _aux_block and just searched for the state, then control, then
            # auxiliary name we would have gotten rk(1)(1) instead of rk(1)
            @test have.args[3] == :(rk__1_ = alpha_.*y__1_./k__1_)
            @test have.args[4] == :(w__1_ = (1.-alpha_).*y__1_./n__1_)
        end

        @testset "negative shift" begin
            have = Dolo._aux_block(sm, -2)
            @test have.head == :block
            @test have.args[1] == :(y_m2_ = z_m2_.*k_m2_.^alpha_.*n_m2_.^(1.-alpha_))
            @test have.args[2] == :(c_m2_ = y_m2_ .- i_m2_)
            @test have.args[3] == :(rk_m2_ = alpha_.*y_m2_./k_m2_)
            @test have.args[4] == :(w_m2_ = (1.-alpha_).*y_m2_./n_m2_)
        end
    end

    @testset "compiled equations" begin
        # the function gets compiled, but can't be called
        eval(Dolo, Dolo.compile_equation(sm, :value))

        @test_throws ErrorException value(m)

        # the transition equations were given in the wrong order. Check that
        # we throw an error here
        @test_throws ErrorException Dolo.compile_equation(sm, :transition)

        # we specified one expecation variable, but two expectation eqns
        @test_throws ErrorException Dolo.compile_equation(sm, :expectation)

        eval(Dolo, Dolo.compile_equation(sm, :auxiliary))

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

        @test_throws MethodError auxiliary(m)
        @test_throws MethodError auxiliary(m, s)
        @test_throws MethodError auxiliary(m, s, x)

        # test allocating version
        @test @inferred(auxiliary(m, s, x, p)) == want

        # test non-allocating version
        @inferred auxiliary!(out, m, s, x, p)
        @test out == want

        # test vectorized version
        out_mat = zeros(5, 4)
        want_mat = [want want want want want]'
        @test @inferred(auxiliary(m, [s s s s s]', x, p)) == want_mat
        @test @inferred(auxiliary(m, s, [x x x x x]', p)) == want_mat

        # non-allocating vectorized version
        @inferred auxiliary!(out_mat, m, s, [x x x x x]', p)
        @test out_mat == want_mat

        # non-allocating vectorized version
        @inferred auxiliary!(out_mat, m, [s s s s s]', x, p)
        @test out_mat == want_mat

        # test errors for wrong size of out
        @test_throws DimensionMismatch auxiliary!(zeros(3), m, s, x, p)
        @test_throws DimensionMismatch auxiliary!(zeros(5), m, s, x, p)
        @test_throws DimensionMismatch auxiliary!(zeros(2, 4), m, [s s s]', x, p)
    end

end
