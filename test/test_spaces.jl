using Test

@testset verbose=true "Cartesian Spaces" begin

@testset verbose=true "Discretization" begin

for d=1:3 

    @testset "Dimension $d" begin

        NN = 20

        vars =  tuple( (Symbol("d$i") for i in 1:d)... )
        pairs = (vars[i]=>[0.0, 1.0*i] for i in 1:d)
        dvars = Dict( pairs )
        cs = Dolo.CSpace(;
            pairs...
        )

        @test isbits(cs)

        @test Dolo.variables(cs) == vars


        # discretization using default values

        cg = Dolo.discretize(cs)
        @test Dolo.size(cg) == tuple( (20 for i=1:d)... )


        # discretization using special values

        dims = tuple( ( 10+i for i=1:d)...  )

        dis = Dict( vars[i]=>dims[i] for i=1:d)
        cg = Dolo.discretize(cs; dis...)
        @test Dolo.size(cg) == dims

        # check it is order insensitive

        dis = Dict( vars[i]=>dims[i] for i=d:-1:1)
        cg = Dolo.discretize(cs; dis...)
        @test Dolo.size(cg) == dims




# f(cg,v) = begin e = sum( sum( e for e in cg ) ); v[1] = e ; nothing end
# res = zeros(1)

# f(cg, res)
# @time f(cg, res)

end

end


end

end