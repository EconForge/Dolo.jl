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

    @testset "_to_expr" begin
        @test Dolo._to_expr("foo") == Expr(:block, :foo)
        @test Dolo._to_expr(100) == Expr(:block, 100)
        @test Dolo._to_expr(:bar) == Expr(:block, :bar)
        @test Dolo._to_expr(:(x+y)) == :(x+y)
    end

    @testset "_expr_or_number" begin
        @test Dolo._expr_or_number(100) == 100
        @test Dolo._expr_or_number(:(x+1)) == :(x+1)
        @test Dolo._expr_or_number(:x) == Expr(:block, :x)
        @test Dolo._expr_or_number("foo") == Expr(:block, :foo)
    end

    @testset "inf_to_Inf" begin
        # numbers and symbols other than :inf go right through
        @test Dolo.inf_to_Inf(42.42) == 42.42
        @test Dolo.inf_to_Inf(:k) == :k

        # :inf and :(-inf) go to Inf and -Inf, respectively
        @test Dolo.inf_to_Inf(:inf) == Inf
        @test Dolo.inf_to_Inf(:(-inf)) == -Inf

        # other expressions are unchanged
        @test Dolo.inf_to_Inf(:(-z)) == :(-z)
    end

    @testset "solve_triangular_system" begin
        # very simple case
        d1 = Dict(:x => 1.0, :y => :(x+1))
        sol1 = Dolo.OrderedDict{Symbol,Number}(:x => 1.0, :y => 2.0)
        @test Dolo.solve_triangular_system(d1) == sol1

        # fully specified numerical dict
        d2 = Dict(:x => 1.0, :y => 2.0)
        sol2 = Dolo.OrderedDict{Symbol,Number}(:x => 1.0, :y => 2.0)
        @test Dolo.solve_triangular_system(d2) == sol2

        # unknown variable w
        d3 = Dict(:x => 1.0, :y => :(x+1), :z=> :(100*w))
        @test_throws ErrorException Dolo.solve_triangular_system(d3)

        # non-triangular system
        d4 = Dict(:x => :(y+1), :y => :(x-1))
        @test_throws ErrorException Dolo.solve_triangular_system(d4)


    end

end
