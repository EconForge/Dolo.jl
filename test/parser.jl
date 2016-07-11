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
            @test Dolo._parse(s) == Symbol(string(s), "_")
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
    calibration
    definitions
end

MockSymbolic(;syms=nothing, eq=nothing, c=nothing, defs=nothing) =
    MockSymbolic{:foobar,:dtcscc}(syms, eq, c, defs)

function Dolo.Options(sm::MockSymbolic, calib::Dolo.ModelCalibration)
    Dolo.Options(Dolo.Cartesian([1.0], [1.0], [1]), Dolo.Normal(eye(2)),
                 Dict{Symbol,Any}())
end

function Dolo.DTCSCCModel{id}(sm::MockSymbolic{id})
    calib = ModelCalibration(FlatCalibration(Dolo.OrderedDict()),
                             GroupedCalibration(Dict()),
                             Dict(),
                             OrderedDict())
    opts = Dolo.Options(sm, calib)
    DTCSCCModel{id}(sm, calib, opts, :dtcscc, "foobar", "boo.yaml")
end

@testset "Applying definitions" begin

    sm = MockSymbolic(;syms=Dict(:states => [:z, :k, :n, :i]),
                       eq=Dict(:transition => [:(y + z*(rk - w)), :(y/z)]),
                       c=Dict(:a=>1.0, :b=>2.0),
                       defs=Dict(:y => parse("z*k^alpha*n^(1-alpha)"),
                                 :c => parse("y - i"),
                                 :rk => parse("alpha*y/k"),
                                 :w => parse("(1-alpha)*y/n")))

    @testset "_is_time_shift" begin
        @test !Dolo._is_time_shift(sm, :(sin(1)))
        @test !Dolo._is_time_shift(sm, :(x(2)))
        @test !Dolo._is_time_shift(sm, :(a(b)))
        @test !Dolo._is_time_shift(sm, :(a(1.0)))
        @test Dolo._is_time_shift(sm, :(a(1)))

        @test Dolo._is_time_shift(:(sin(1)))
        @test Dolo._is_time_shift(:(x(2)))
        @test !Dolo._is_time_shift(:(a(b)))
        @test !Dolo._is_time_shift(:(a(1.0)))
        @test Dolo._is_time_shift(:(a(1)))
    end

    @testset "_replace_with_shift" begin
        rws = Dolo._replace_with_shift

        @test rws(sm, 1, :c, 1) == 1
        @test rws(sm, 1, :c, 1432143) == 1
        @test rws(sm, :c, :c, 1) == :(c(1))
        @test rws(sm, :c, :c, -1) == :(c(-1))
        @test rws(sm, :x, :c, 1) == :x
        @test rws(sm, :(a*y - a), :a, 1) == :(a(1)*y - a(1))
        @test rws(sm, :(a*y - a), :y, -1) == :(a*y(-1) - a)
        @test rws(sm, :(a*y - a), :b, 10) == :(a*y - a)
    end

    @testset "_call_definition" begin
        _cd = Dolo._call_definition

        want = :(z(0) * k(0) ^ alpha * n(0) ^ (1 - alpha))
        @test _cd(sm, sm.definitions[:y], 0) == want

        want = :(z(-1) * k(-1) ^ alpha * n(-1) ^ (1 - alpha))
        @test _cd(sm, sm.definitions[:y], -1) == want

        # test recursive application
        want = :(z(0) * k(0) ^ alpha * n(0) ^ (1 - alpha) - i(0))
        @test _cd(sm, sm.definitions[:c], 0) == want

        want = :(z(-1) * k(-1) ^ alpha * n(-1) ^ (1 - alpha) - i(-1))
        @test _cd(sm, sm.definitions[:c], -1) == want
    end

    @testset "_apply_definitions" begin
        _ad = Dolo._apply_definitions

        want1 = :(z(0) * k(0) ^ alpha * n(0) ^ (1 - alpha) + z * ((alpha * (z(0) * k(0) ^ alpha * n(0) ^ (1 - alpha))) / k(0) - ((1 - alpha) * (z(0) * k(0) ^ alpha * n(0) ^ (1 - alpha))) / n(0)))
        want2 = :((z(0) * k(0) ^ alpha * n(0) ^ (1 - alpha)) / z)

        @test _ad(sm, sm.equations[:transition][1]) == want1
        @test _ad(sm, sm.equations[:transition][2]) == want2
        @test _ad(sm, sm.equations[:transition]) == [want1, want2]
    end
end

@testset "compiler" begin

    @testset "_param_block" begin
        sm = MockSymbolic(;syms=Dict(:parameters => [:a, :b, :foobar]))
        have = Dolo._param_block(sm)
        @test have.head == :block
        @test have.args[1] == :(a_ = _unpack_var(p, 1))
        @test have.args[2] == :(b_ = _unpack_var(p, 2))
        @test have.args[3] == :(foobar_ = _unpack_var(p, 3))
    end

    @testset "_single_arg_block" begin
        sm = MockSymbolic(;syms=Dict(:states => [:z, :k]))

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
        sm = MockSymbolic(;syms=Dict())
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
    defs = Dict(:y => parse("z*k^alpha*n^(1-alpha)"),
                :c => parse("y - i"),
                :rk => parse("alpha*y/k"),
                :w => parse("(1-alpha)*y/n"))
    sm = MockSymbolic(;syms=Dict(:states => [:z, :k],
                                 :shocks => [:ɛz],
                                 :controls => [:i, :n],
                                 :expectations => [:Ez],
                                 :parameters => [:fizz, :buzz, :alpha]),
                      eq=Dict(:value => Expr[],
                              :arbitrage => [:(1 - fizz*(c/c(1))^(buzz)*(1-alpha+rk(1))),
                                             :(z*n^fizz*c^buzz - w)],
                              :transition => [:(k = k(-1) + 1), :(z = z(-1))],
                              :expectation => [:(Ez = z(1)), :(Ek = k(1))]),
                      defs=defs)

    m = DTCSCCModel(sm)

    @testset "compiled equations" begin
        # the function gets compiled, but can't be called
        eval(Dolo, Dolo.compile_equation(sm, :value))

        @test_throws ErrorException value(m)

        # the transition equations were given in the wrong order. Check that
        # we throw an error here
        @test_throws ErrorException Dolo.compile_equation(sm, :transition)

        # we specified one expecation variable, but two expectation eqns
        @test_throws ErrorException Dolo.compile_equation(sm, :expectation)

        # now see if the functions were compiled properly
        eval(Dolo, Dolo.compile_equation(sm, :arbitrage))

        z = 1.0
        k = 0.25
        i = 100.0
        n = 10.0
        alpha = 0.33
        fizz = 1.1
        buzz = 2
        s = [z, k]
        x = [i, n]
        E = [0.0]
        S = s
        X = x
        p = [fizz, buzz, alpha]

        y = z*k^alpha*n^(1-alpha)
        c = y-i
        rk = alpha*y/k
        w = (1-alpha)*y/n

        want = [1 - fizz*(c/c)^(buzz)*(1-alpha+rk),
                z*n^fizz*c^buzz - w]
        out = zeros(2)

        @test_throws MethodError arbitrage(m)
        @test_throws MethodError arbitrage(m, s)
        @test_throws MethodError arbitrage(m, s, x)
        @test_throws MethodError arbitrage(m, s, x, E)
        @test_throws MethodError arbitrage(m, s, x, E, S)
        @test_throws MethodError arbitrage(m, s, x, E, S, X)

        # test allocating version
        @test @inferred(arbitrage(m, s, x, E, S, X, p)) == want

        # test non-allocating version
        @inferred arbitrage!(out, m, s, x, E, S, X, p)
        @test out == want

        # test vectorized version
        out_mat = zeros(5, 2)
        want_mat = [want want want want want]'
        @test @inferred(arbitrage(m, [s s s s s]', x, E, S, X, p)) == want_mat
        @test @inferred(arbitrage(m, s, [x x x x x]', E, S, X, p)) == want_mat

        # non-allocating vectorized version
        @inferred arbitrage!(out_mat, m, s, [x x x x x]', E, S, X, p)
        @test out_mat == want_mat

        # non-allocating vectorized version
        @inferred arbitrage!(out_mat, m, [s s s s s]', x, E, S, X, p)
        @test out_mat == want_mat

        # test errors for wrong size of out
        @test_throws DimensionMismatch arbitrage!(zeros(3), m, s, x, E, S, X, p)
        @test_throws DimensionMismatch arbitrage!(zeros(5), m, s, x, E, S, X, p)
        @test_throws DimensionMismatch arbitrage!(zeros(2, 4), m, [s s s]', x, E, S, X, p)
    end

end
