@testset "_handle_arbitrage" begin
    arb = ["x + y + z + w",
           "x + y + z + w | y <= 100",
           "x + y + z + w | 0 <= z",
           "x + y + z + w | -42.0 <= w <= 42.0",
           ]
    controls = [:x, :y, :z, :w]

    lb, ub, a = Dolo._handle_arbitrage(arb, controls)

    @test length(lb) == 4
    @test length(ub) == 4
    @test length(a) == 4

    @test lb[1] == Expr(:block, -Inf)
    @test lb[2] == Expr(:block, -Inf)
    @test lb[3] == Expr(:block, 0)
    @test lb[4] == Expr(:block, -42.0)

    @test ub[1] == Expr(:block, Inf)
    @test ub[2] == Expr(:block, 100)
    @test ub[3] == Expr(:block, Inf)
    @test ub[4] == Expr(:block, 42.0)


    ex = :(x + y + z + w)
    for i in 1:4
        @test a[i] == ex
    end

    # now test that errors are thrown properly

    @testset "error if control in control is wrong in one sided comp." begin
        @test_throws ErrorException Dolo._handle_arbitrage(["x+y|i<=10"], [:y])
    end

    @testset "error if control in control is wrong in two sided comp." begin
        eqn = ["x+y|0 <= i<= 10"]
        @test_throws ErrorException Dolo._handle_arbitrage(["x+y|i<=10"], [:y])
    end

    @testset "error if multiple |" begin
        eqn = ["x+y|i<=10 | 0 <= i"]
        @test_throws ErrorException Dolo._handle_arbitrage(eqn, [:i])
    end

    @testset "error if 3 <=" begin
        eqn = ["x+y| 0 <= 1 <= i <= 100"]
        @test_throws ErrorException Dolo._handle_arbitrage(eqn, [:i])
    end

    @testset "error if zero <= (maybe some other inequality was used)" begin
        eqn = ["x+y| 0 < y < 100"]
        @test_throws ErrorException Dolo._handle_arbitrage(eqn, [:i])
    end
end

@testset "_parse" begin
    @testset "numbers" begin
        for T in (Float16, Float32, Float64, Int8, Int16, Int32, Int64)
            x = rand(T)
            @test _parse(x) == x
        end
    end

    @testset "symbols" begin
        for i=1:10
            s = gensym()
            @test _parse(s) == symbol(string(s), "_")
        end
    end

    @testset "x_(shift_Integer)" begin
        for i=1:10, T in (Int8, Int16, Int32, Int64)
            @test _parse(string("x(", T(i), ")")) == symbol("x__$(i)_")
            @test _parse(string("x(", T(-i), ")")) == symbol("x_m$(i)_")
        end

        # only add underscore to naems when shift is 0
        @test _parse("x(0)") == :x_
    end

    @testset "other function calls" begin
        @testset "one argument" begin
            @test _parse("sin(x)") == :(sin(x_))
            @test _parse("sin(x(-1))") == :(sin(x_m1_))
            @test _parse("foobar(x(2))") == :(foobar(x__2_))
        end

        @testset "two arguments" begin
            @test _parse("dot(x, y(1))") == :(dot(x_, y__1_))
            @test _parse("plot(x(-1), y)") == :(plot(x_m1_, y_))
            @test _parse("bingbong(x(2), y)") == :(bingbong(x__2_, y_))
        end

        @testset "more args" begin
            for i=3:10
                ex = Expr(:call, :(+), [:(x($j)) for j in 1:i]...)
                want = Expr(:call, :(+), [symbol("x__", j, "_") for j in 1:i]...)
                @test _parse(ex) == want
            end
        end

        @testset "throws errors when unsupported" begin
            @test_throws ErrorException _parse("x+y | i <= 100")
        end
    end

    @testset "Expr(:(=), ...)" begin
        @testset "without targets" begin
            @test _parse(:(x = y)) == :(y_ - x_)
        end

        @testset "with targets" begin
            @test _parse(:(x = log(y(-1))); targets=[:x]) == :(x_ = log(y_m1_))
            @test_throws ErrorException _parse(:(x = y); targets=[:y])
        end
    end
end
