

@testset "testing model_types" begin
    #=
    I got this by doing

    ```
    src = load_file("/Users/sglyon/src/Python/dolo/examples/models/rbc.yaml")
    repr(src)
    ```
    And then doing some formatting
    =#

    rbc_dict = Dict{Any,Any}(
        "symbols"=>Dict{Any,Any}(
            "parameters"=>Any["beta","sigma","eta","chi","delta","alpha","rho","zbar","sig_z"],
            "controls"=>Any["i","n"],
            "values"=>Any["V"],
            "auxiliaries"=>Any["y","c","rk","w"],
            "states"=>Any["z","k"],
            "shocks"=>Any["e_z"]
        ),
        "name"=>"Real Business Cycle",
        "calibration"=>Dict{Any,Any}(
            "c"=>"y - i","zbar"=>1,
            "V"=>"log(c)/(1-beta)",
            "delta"=>0.025,"sigma"=>1,"chi"=>"w/c^sigma/n^eta",
            "phi"=>1,"z"=>"zbar",
            "rk"=>"1/beta-1+delta",
            "i"=>"delta*k",
            "y"=>"z*k^alpha*n^(1-alpha)",
            "sig_z"=>0.016,"w"=>"(1-alpha)*z*(k/n)^(alpha)",
            "alpha"=>0.33,"k"=>"n/(rk/alpha)^(1/(1-alpha))",
            "eta"=>1,
            "rho"=>0.8,
            "beta"=>0.99,
            "n"=>0.33
        ),
        "options"=>Dict{Any,Any}(
            "Approximation"=>Dict{Any,Any}(
                "orders"=>Any[10,50],
                "b"=>Any["1+2*sig_z","k*1.1"],
                "a"=>Any["1-2*sig_z","k*0.9"]
                )
            ),
        "equations"=>Dict{Any,Any}(
            "arbitrage"=>Any[
                "1 - beta*(c/c(1))^(sigma)*(1-delta+rk(1))   | 0 <= i <= inf",
                "chi*n^eta*c^sigma - w                       | 0 <= n <= inf"
            ],
            "auxiliary"=>Any[
                "y = z*k^alpha*n^(1-alpha)","c = y - i",
                "rk = alpha*y/k",
                "w = (1-alpha)*y/n"
            ],
            "value"=>Any["V = log(c) + beta*V(1)"],
            "transition"=>Any[
                "z = (1-rho)*zbar + rho*z(-1) + e_z",
                "k = (1-delta)*k(-1) + i(-1)"
            ]
        ),
        "distribution"=>Dict{Any,Any}(
            "Normal"=>Any[Any["sig_z**2"]]
        )
    )

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

    @testset "SymbolicModel" begin
        sm = SymbolicModel(rbc_dict, :dtcscc, "rbc.yaml")

        @testset "constructor" begin

            # construct the object, and now test to make sure that things are how
            # they should be

            @test length(keys(sm.symbols)) == 6
            np, nc, nv, na, ns, nz = 9, 2, 1, 4, 2, 1
            @test length(sm.symbols[:parameters]) == np
            @test length(sm.symbols[:controls]) == nc
            @test length(sm.symbols[:values]) == nv
            @test length(sm.symbols[:auxiliaries]) == na
            @test length(sm.symbols[:states]) == ns
            @test length(sm.symbols[:shocks]) == nz

            # test order of keys
            @test collect(keys(sm.symbols)) == map(symbol, Dolo.RECIPES[:dtcscc][:symbols])

            # test that we got the right number of equations for each group
            @test length(sm.equations[:arbitrage]) == nc
            @test length(sm.equations[:value]) == nv
            @test length(sm.equations[:transition]) == ns
            @test length(sm.equations[:auxiliary]) == na
            @test length(sm.equations[:controls_lb]) == nc
            @test length(sm.equations[:controls_ub]) == nc

            name_order = vcat([sm.symbols[symbol(k)]
                              for k in Dolo.RECIPES[:dtcscc][:symbols]]...)
            @test collect(keys(sm.calibration)) == name_order
            @test sm.options == Dolo._symbol_dict(rbc_dict["options"])
            @test sm.distribution == Dolo._symbol_dict(rbc_dict["distribution"])
        end

        @test model_type(sm) == sm.model_type == :dtcscc
        @test name(sm) == sm.name == "Real Business Cycle"
        @test filename(sm) == sm.filename == "rbc.yaml"
    end

    @testset "NumericModel" begin
        sm = SymbolicModel(rbc_dict, :dtcscc, "rbc.yaml")
        m = DTCSCCModel(sm)

        # check that options and distribution were "numericized" properly
        @test all([isa(x, Number) for x in m.options[:Approximation][:orders]])
        @test all([isa(x, Number) for x in m.options[:Approximation][:b]])
        @test all([isa(x, Number) for x in m.options[:Approximation][:a]])
        @test all([isa(x, Number) for x in m.distribution[:Normal]])
        @test isa(m.distribution[:Normal], Matrix)

        @test model_type(m) == m.model_type == :dtcscc
        @test name(m) == m.name == "Real Business Cycle"
        @test filename(m) == m.filename == "rbc.yaml"
    end

    @testset "DTCSCCfunctions" begin
        sm = SymbolicModel(rbc_dict, :dtcscc, "rbc.yaml")
        m = DTCSCCModel(sm)

        # we passed in the steady state, so just make sure compiled functions
        # work out that way
        mc = m.calibration
        p, s0, x0, a0, v0 = mc[:parameters, :states, :controls, :auxiliaries, :values]
        e_ = zeros(mc[:shocks])
        ns = length(s0)
        nx = length(x0)

        funcs = m.functions

        @test maxabs(zeros(x0) - evaluate(funcs.arbitrage, s0, x0, e_, s0, x0, p)) < 1e-13
        @test maxabs(s0 - evaluate(funcs.transition, s0, x0, e_, p)) < 1e-13
        @test maxabs(a0 - evaluate(funcs.auxiliary, s0, x0, p)) < 1e-13
        @test maxabs(v0 - evaluate(funcs.value, s0, x0, s0, x0, v0, p)) < 1e-13
        @test maxabs([0.0, 0.0] - evaluate(funcs.controls_lb, s0, p)) < 1e-13
        @test [Inf, Inf] == evaluate(funcs.controls_ub, s0, p)

        # these two aren't implemented for the model above
        @test_throws ErrorException evaluate(funcs.direct_response, s0, [0.0], p)
        @test_throws ErrorException evaluate(funcs.expectation, s0, x0, p)

        # Test mutating versions
        res = ones(x0)
        evaluate!(funcs.arbitrage, s0, x0, e_, s0, x0, p, res)
        @test maxabs(zeros(x0) - res) < 1e-13

        s1 = ones(s0)
        evaluate!(funcs.transition, s0, x0, e_, p, s1)
        @test maxabs(s0 - s1) < 1e-13

        a1 = ones(a0)
        evaluate!(funcs.auxiliary, s0, x0, p, a1)
        @test maxabs(a0 - a1) < 1e-13

        v1 = ones(v0)
        evaluate!(funcs.value, s0, x0, s0, x0, v0, p, v1)
        @test maxabs(v0 - v1) < 1e-13

        bounds = ones(x0)
        evaluate!(funcs.controls_ub, s0, p, bounds)
        @test bounds == [Inf, Inf]

        evaluate!(funcs.controls_lb, s0, p, bounds)
        @test bounds == [0.0, 0.0]
    end

end
