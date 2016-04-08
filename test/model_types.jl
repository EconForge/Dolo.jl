@testset "testing model_types" begin

    @testset "ModelCalibration" begin
        function new_mc()
            flat = FlatCalibration(:k=>8.5, :z=>0.5, :i=>1.1)
            grouped = GroupedCalibration(:states=>[8.5, 0.5],
                                                   :controls=>[1.1])
            symbol_table = OrderedDict{Symbol,Tuple{Symbol,Int}}()
            symbol_table[:k] = (:states, 1)
            symbol_table[:z] = (:states, 2)
            symbol_table[:i] = (:controls, 1)
            symbol_groups = OrderedDict{Symbol,Vector{Symbol}}()
            symbol_groups[:states] = [:k, :z]
            symbol_groups[:controls] = [:i]
            mc = ModelCalibration(flat, grouped, symbol_table, symbol_groups)
        end
        mc = new_mc()
        mc2 = copy(mc)  # shallow copy. underlying arrays are same
        mc3 = deepcopy(mc)  # deep copy. underlying arrays are different

        @testset "FlatCalibration" begin
            d = OrderedDict(:x=>1.0, :y=>2.0)

            # test constructors
            fc1 = FlatCalibration(d)
            fc2 = FlatCalibration(:x=>1.0, :y=>2.0)
            @test fc1 == fc2
            @test !(fc1 === fc2)

            # getindex
            @test 1.0 == @inferred getindex(fc1, :x)
            @test [1.0, 2.0] == @inferred getindex(fc1, :x, :y)

            # setindex
            fc1[:z] = 1.0
            @test fc1[:z] == 1.0

            # conversion to Float64
            for T in (Float16, Float32, Float64, Int8, Int16, Int32, Int64,
                      Rational{Int}, BigFloat, BigInt)
                fc1[:q] = one(T)
                @test Float64(1.0) == @inferred getindex(fc1, :q)
            end

            # multiple values from tuple
            fc1[:a, :b, :c] = 4, 5, 6
            @test 4.0 == @inferred getindex(fc1, :a)
            @test 5.0 == @inferred getindex(fc1, :b)
            @test 6.0 == @inferred getindex(fc1, :c)
            @test_throws DimensionMismatch setindex!(fc1, (4, 5), :a, :b, :c)

            # multiple values from (poorly typed) Vector
            fc1[:a, :b, :c] = Any[7, 8, 9]
            @test 7.0 == @inferred getindex(fc1, :a)
            @test 8.0 == @inferred getindex(fc1, :b)
            @test 9.0 == @inferred getindex(fc1, :c)
            @test_throws DimensionMismatch setindex!(fc1, [7, 8], :a, :b, :c)
        end

        @testset "GroupedCalibration" begin
            d = Dict(:x=>[1.0, 2.0], :y=>[2.0, 3.0])

            # test constructors
            gc1 = GroupedCalibration(d)
            gc2 = GroupedCalibration(:x=>[1.0, 2.0], :y=>[2.0, 3.0])
            @test gc1 == gc2
            @test !(gc1 === gc2)

            # getindex
            @test [1.0, 2.0] == @inferred getindex(gc1, :x)
            @test Vector{Float64}[[1.0, 2.0], [2.0, 3.0]] == @inferred getindex(gc1, :x, :y)

            # setindex
            gc1[:z] = [1.0]
            @test gc1[:z] == [1.0]

            # existing x has length 2, this is only length 1
            @test_throws DimensionMismatch setindex!(gc1, [1.0], :x)

            # conversion to Float64
            for T in (Float16, Float32, Float64, Int8, Int16, Int32, Int64,
                      Rational{Int}, BigFloat, BigInt)
                gc1[:q] = [one(T)]
                @test [Float64(1.0)] == @inferred getindex(gc1, :q)
            end

            # multiple values from tuple
            gc1[:a, :b, :c] = [4], [5], [6]
            @test [4.0] == @inferred getindex(gc1, :a)
            @test [5.0] == @inferred getindex(gc1, :b)
            @test [6.0] == @inferred getindex(gc1, :c)
            @test_throws DimensionMismatch setindex!(gc1, ([4], [5]), :a, :b, :c)

        end

        @testset "copy/deepcopy" begin

            # triple = effectively checks if two arrays point to same memory
            # address
            @test mc2.grouped[:states] === mc.grouped[:states]
            @test !(mc3.grouped[:states] === mc.grouped[:states])
        end

        @testset "getindex/setindex!" begin
            # getindex
            @test [8.5, 0.5] == @inferred getindex(mc, :states)
            @test Vector{Float64}[[8.5, 0.5], [1.1]] == @inferred getindex(mc, :states, :controls)

            # test setindex! for a group
            mc3[:states] = [42.0, 43.0]
            @test [42.0, 43.0] == @inferred getindex(mc3, :states)
            @test 42.0 == mc3.flat[:k]
            @test 43.0 == mc3.flat[:z]
            @test_throws DimensionMismatch setindex!(mc3, rand(3), :states)

            # test setindex! for multiple groups
            mc3[:states, :controls] = [43.0, 42.0], [1.5]
            @test [43.0, 42.0] == @inferred getindex(mc3, :states)
            @test [1.5] == @inferred getindex(mc3, :controls)
            @test 43.0 == mc3.flat[:k]
            @test 42.0 == mc3.flat[:z]
            @test 1.5 == mc3.flat[:i]
            @test_throws DimensionMismatch setindex!(mc3, (rand(2),), :states, :controls)
            @test_throws DimensionMismatch setindex!(mc3, (rand(3),rand(1)), :states, :controls)
        end

        @testset "eval_with" begin
            # test _replace_me
            @test Dolo._replace_me(mc, :(+)) == :(+)
            @test Dolo._replace_me(mc, :(bing)) == :(bing)
            @test Dolo._replace_me(mc, 100) == 100
            @test Dolo._replace_me(mc, :k) == 8.5

            @test eval_with(mc, "k+1") == 9.5
            @test eval_with(mc, "k+i") == 9.6

            # test expression version
            @test eval_with(mc, :(k+i)) == 9.6

            # test that we can define temporaries that aren't created in this
            # scope
            @test 9.6 == eval_with(mc, quote
                foobar = k
                foobar + i
                end)
            @test !isdefined(current_module(), :foobar)

        end
    end
end
