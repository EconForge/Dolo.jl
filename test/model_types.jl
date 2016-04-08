@testset "testing model_types" begin

    @testset "ModelCalibration" begin
        function new_mc()
            flat = OrderedDict{Symbol,Float64}(:k=>8.5, :z=>0.5, :i=>1.1)
            grouped = Dict{Symbol,Vector{Float64}}(:states=>[8.5, 0.5],
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

        @testset "copy/deepcopy" begin

            # triple = effectively checks if two arrays point to same memory
            # address
            @test mc2.grouped[:states] === mc.grouped[:states]
            @test !(mc3.grouped[:states] === mc.grouped[:states])
        end

        @testset "getindex/setindex!" begin
            # getindex
            @test 8.5 == mc[:k] == @inferred getindex(mc, :k)
            @test 8.5 == mc.flat[:k] == @inferred getindex(mc.flat, :k)
            @test [8.5, 0.5] == mc["states"] == @inferred getindex(mc, "states")
            @test [8.5, 0.5] == mc.grouped[:states] == @inferred getindex(mc.grouped, :states)

            # setindex!. Do this on mc3 so as not to change mc. Verify that
            # arrays in mc are unchanged after this
            mc3[:k] = 5.0

            @test mc3[:k] == 5.0
            @test mc3["states"] == [5.0, 0.5]
            @test mc[:k] == 8.5
            @test mc["states"] == [8.5, 0.5]

            # try same with mc2 and show that mc[:k] is unchanged,
            # but mc["states"] is
            mc2[:k] = 5.0
            @test mc2[:k] == 5.0
            @test mc2["states"] == [5.0, 0.5]
            @test mc[:k] == 8.5
            @test mc["states"] == [5.0, 0.5]

            # make sure we throw if setindex! doesn't have matching sizes
            @test_throws ErrorException setindex!(mc3, rand(4), :k, :i, :z)

            # now get back our original mc
            mc = new_mc()

            # test vector version
            @test [8.5, 1.1] == @inferred getindex(mc, :k, :i)
            @test Vector{Float64}[[8.5, 0.5], [1.1]] == @inferred getindex(mc, "states", "controls")

            # test setindex! for a group
            mc3["states"] = [42.0, 43.0]
            @test [42.0, 43.0] == @inferred getindex(mc3, "states")
            @test 42.0 == mc3[:k]
            @test 43.0 == mc3[:z]
            @test_throws ErrorException setindex!(mc3, rand(3), "states")
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
