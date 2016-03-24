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
