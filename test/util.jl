@testset "testing util" begin

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

    @testset "_expr_or_number" begin
        @test Dolo._expr_or_number(100) == 100
        @test Dolo._expr_or_number(:(x+1)) == :(x+1)
        @test Dolo._expr_or_number(:x) == Expr(:block, :x)
        @test Dolo._expr_or_number("foo") == Expr(:block, :foo)
    end    
end
