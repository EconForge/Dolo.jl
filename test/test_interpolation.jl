using Test
using StaticArrays

@testset verbose=true "Cartesian Interpolation" begin

@testset verbose=true "Interpolation" begin

    for d in 1:3 

        @testset verbose=true "Dimension $d" begin

            NN = 20

            vars =  tuple( (Symbol("d$i") for i in 1:d)... )
            pairs = (vars[i]=>(0.0, 1.0*i) for i in 1:d)
            dvars = Dict( pairs )
            cs = Dolo.CSpace(;
                pairs...
            )

            sg = Dolo.discretize(cs)

            coeffs = [(i+1)^2 for i=1:d]
            coeffs_2 = [(i+1)^3 for i=1:d]

            fl(x) = SVector( sum(x .* coeffs), sum(x.* coeffs_2) )
            fc(x) = SVector( sum(x .* coeffs + x.^2 .* coeffs), sum(x.* coeffs_2  + x.^3 .* coeffs) )
            

            function testfun(dfun, sg)
                if sum(sum(dfun(e) for e in sg)) <0.0
                    println("Oups")
                else
                    nothing
                end
            end

            @testset verbose=true "Linear" begin

                values = [fl(e) for e in sg]
                gvec = Dolo.GVector(sg, values)

                dfun = Dolo.DFun(cs, gvec)

                _values = [dfun(e) for e in sg]

                d = _values - values
                @test isapprox(values, _values)
                testfun(dfun, sg)
                @test 0== @allocated testfun(dfun, sg)

            end

            @testset verbose=true "Cubic" begin

                values = [fc(e) for e in sg]
                gvec = Dolo.GVector(sg, values)

                dfun = Dolo.DFun(cs, gvec; interp_mode=:cubic)

                _values = [dfun(e) for e in sg]

                d = _values - values
                @test isapprox(values, _values)
                
                testfun(dfun, sg)
                @test 0== @allocated testfun(dfun, sg)
                
            end

        end

    end

end


end
