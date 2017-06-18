
@testset "testing model_types" begin
    #=
    I got this by doing

    ```
    using YAML: load_file
    using Dolo
    funcs = Dict("!Cartesian" => (c, n) -> Dolo.construct_type_map(:Cartesian, c, n),
             "!Normal" => (c, n) -> Dolo.construct_type_map(:Normal, c, n))
    src = load_file("/Users/sglyon/src/Python/dolo/examples/models/rbc.yaml",
                    funcs)
    print(repr(src))
    ```
    And then doing some formatting
    =#
    #
    # rbc_dict = Dict{Any,Any}(
    #     "symbols"=>Dict{Any,Any}(
    #         "parameters"=>Any["beta","sigma","eta","chi","delta","alpha","rho","zbar","sig_z"],
    #         "controls"=>Any["i","n"],
    #         "values"=>Any["V"],
    #         "states"=>Any["k"],
    #         "exogenous"=>Any["z"]
    #     ),
    #     "model_type"=>"dtcc",
    #     "name"=>"Real Business Cycle",
    #     "calibration"=>Dict{Any,Any}(
    #         "c"=>"y - i",
    #         "zbar"=>1,
    #         "V"=>"log(c)/(1-beta)",
    #         "c_i"=>1.5,
    #         "delta"=>0.025,
    #         "sigma"=>1,
    #         "chi"=>"w/c^sigma/n^eta",
    #         "phi"=>1,
    #         "z"=>"zbar",
    #         "rk"=>"1/beta-1+delta",
    #         "i"=>"delta*k",
    #         "y"=>"z*k^alpha*n^(1-alpha)",
    #         "sig_z"=>0.016,
    #         "w"=>"(1-alpha)*z*(k/n)^(alpha)",
    #         "alpha"=>0.33,
    #         "k"=>"n/(rk/alpha)^(1/(1-alpha))",
    #         "eta"=>1,
    #         "rho"=>0.8,
    #         "beta"=>0.99,
    #         "c_y"=>0.5,
    #         "n"=>0.33),
    #
    #     "options"=>Dict{Any,Any}(
    #         # "distribution"=>Dict{Symbol,Any}(
    #         "exogenous"=>Dict{Any,Any}(
    #             :rho => 0.9,
    #             :sigma => "sig_z",
    #             :tag => :AR1
    #           ),
    #         #     :sigma=>Any[Any["sig_z**2"]],
    #         #     :tag=>:Normal),
    #         "grid"=>Dict{Symbol,Any}(
    #             :orders=>Any[10],
    #             :b=>Any["k*1.1"],
    #             :a=>Any["k*0.9"],
    #             :tag=>:Cartesian
    #           )
    #     ),
    #     "equations"=>Dict{Any,Any}(
    #         "arbitrage"=>Any["1 - beta*(c/c(1))^(sigma)*(1-delta+rk(1))  | 0 <= i <= inf","chi*n^eta*c^sigma - w                      | 0 <= n <= inf"],
    #         "value"=>Any["V = log(c) + beta*V(1)"],
    #         "transition"=>Any["k = (1-delta)*k(-1) + i(-1)"]
    #     ),
    #     "definitions"=>Dict{Any,Any}(
    #         "c"=>"y - i",
    #         "w"=>"(1-alpha)*y/n",
    #         "rk"=>"alpha*y/k",
    #         "y"=>"z*k^alpha*n^(1-alpha)"
    #     )
    # )
    #
    # rbc_dict = Dolo._symbol_dict(rbc_dict)


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
            for T in (Float16, Float32, Float64, Int8, Int16, Int32, Int,
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
            for T in (Float16, Float32, Float64, Int8, Int16, Int32, Int,
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

            # setindex!. Do this on mc3 so as not to change mc. Verify that
            # arrays in mc are unchanged after this
            mc3[:k] = 5.0

            @test mc3.flat[:k] == 5.0
            @test mc3[:states] == [5.0, 0.5]
            @test mc.flat[:k] == 8.5
            @test mc[:states] == [8.5, 0.5]

            # try same with mc2 and show that mc[:k] is unchanged,
            # but mc["states"] is
            mc2[:k] = 5.0
            @test mc2.flat[:k] == 5.0
            @test mc2[:states] == [5.0, 0.5]
            @test mc.flat[:k] == 8.5
            @test mc[:states] == [5.0, 0.5]

            # make sure we throw if setindex! doesn't have matching sizes
            @test_throws DimensionMismatch setindex!(mc3, (1,2,3,4), :k, :i, :z)

            # now get back our original mc
            mc = new_mc()

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

            @test eval_with(mc, ["k+1", "k+i"]) == [9.5, 9.6]
            @test eval_with(mc, 1.0) == 1.0
            @test eval_with(mc, :k) == 8.5
            @test eval_with(mc, Dict("hi"=>"i", "low"=>"-i")) ==
                    Dict{Symbol,Any}(:hi => 1.1, :low => -1.1)

            # make sure it works with any associatives
            @test eval_with(mc, mc.flat) == mc.flat.d
            @test eval_with(mc, mc.grouped) == mc.grouped.d

        end
    end
    #
    # @testset "SymbolicModel" begin
    #     sm = SymbolicModel(rbc_dict, "rbc_dtcc.yaml")
    #
    #     @testset "constructor" begin
    #
    #         # construct the object, and now test to make sure that things are how
    #         # they should be
    #
    #         @test length(keys(sm.symbols)) == 8
    #         np, nc, nv, na, ns, nz = 9, 2, 1, 1, 1, 1
    #         @test length(sm.symbols[:parameters]) == np
    #         @test length(sm.symbols[:controls]) == nc
    #         @test length(sm.symbols[:values]) == nv
    #         @test length(sm.symbols[:states]) == ns
    #         @test length(sm.symbols[:exogenous]) == nz
    #         @test length(sm.symbols[:rewards]) == 0
    #         @test length(sm.symbols[:expectations]) == 0
    #
    #         # test order of keys
    #         @test collect(keys(sm.symbols)) == map(Symbol, Dolo.RECIPES[:dtcc][:symbols])
    #
    #         # test that we got the right number of equations for each group
    #         @test length(sm.equations[:arbitrage]) == nc
    #         @test length(sm.equations[:value]) == nv
    #         @test length(sm.equations[:transition]) == ns
    #         @test length(sm.equations[:controls_lb]) == nc
    #         @test length(sm.equations[:controls_ub]) == nc
    #         @test !haskey(sm.equations, :rewards)
    #         @test !haskey(sm.equations, :expectations)
    #
    #         name_order = vcat([sm.symbols[Symbol(k)]
    #                           for k in Dolo.RECIPES[:dtcc][:symbols]]...,
    #                           collect(keys(sm.definitions)))
    #         @test collect(keys(sm.calibration)) == name_order
    #     end
    #
    #     @test name(sm) == sm.name == "Real Business Cycle"
    #     @test filename(sm) == sm.filename == "rbc_dtcc.yaml"
    # end
    #
    # @testset "NumericModel" begin
    #     sm = SymbolicModel(rbc_dict, "rbc_dtcc.yaml")
    #     m = NumericModel(sm)
    #
    #     @test name(m) == m.name == "Real Business Cycle"
    #     @test filename(m) == m.filename == "rbc_dtcc.yaml"
    #     @test id(sm) == id(m)
    # end
    #
    # @testset "compiled functions" begin
    #     sm = SymbolicModel(rbc_dict, "rbc_dtcc.yaml")
    #     m = NumericModel(sm)
    #
    #     # we passed in the steady state, so just make sure compiled functions
    #     # work out that way
    #     mc = m.calibration
    #     p, s0, x0, v0 , e_= mc[:parameters, :states, :controls, :values, :exogenous]
    #     ns = length(s0)
    #     nx = length(x0)
    #
    #     @test maximum(abs, zeros(x0) - arbitrage(m, e_, s0, x0, e_, s0, x0, p)) < 1e-13
    #     @test maximum(abs, s0 - transition(m, e_, s0, x0, e_, p)) < 1e-13
    #     @test maximum(abs, v0 - value(m, e_, s0, x0, v0, e_, s0, x0, v0, p)) < 1e-13
    #     @test maximum(abs, [0.0, 0.0] - controls_lb(m, e_, s0, p)) < 1e-13
    #     @test [Inf, Inf] == controls_ub(m, e_, s0, p)
    #
    #     # these two aren't implemented for the model above
    #     @test_throws ErrorException direct_response(m, e_, s0, [0.0], p)
    #     @test_throws ErrorException expectation(m, e_, s0, x0, p)
    #
    #     # Test mutating versions
    #     res = ones(x0)
    #     arbitrage!(m, res, e_, s0, x0, e_, s0, x0, p)
    #     @test maximum(abs, zeros(x0) - res) < 1e-13
    #
    #     s1 = ones(s0)
    #     transition!(m, s1, e_, s0, x0, e_, p)
    #     @test maximum(abs, s0 - s1) < 1e-13
    #
    #     v1 = ones(v0)
    #     value!(m, v1, e_, s0, x0, v0, e_, s0, x0, v0, p)
    #     @test maximum(abs, v0 - v1) < 1e-13
    #
    #     bounds = ones(x0)
    #     controls_ub!(m, bounds, e_, s0, p)
    #     @test bounds == [Inf, Inf]
    #
    #     controls_lb!(m, bounds, e_, s0, p)
    #     @test bounds == [0.0, 0.0]
    # end

end
