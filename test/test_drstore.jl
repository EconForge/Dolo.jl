
import Dolo: DRStore
using StaticArrays

@testset "DR Storage" begin

    @testset "floats"

        vals = rand(150)
        dims = [50,60,40]
        ddr = DRStore(vals, dims)

        ddrc = copy(ddr)

        ddr.flat[3] = 19.341
        @test maximum(abs,ddr.flat-vals)==0.0
        @test maximum(abs,ddr.flat-ddrc.flat)!=0.0
        @test maximum(cat(1, ddr.data...) - ddr.flat)==0.0


        dd = [rand(4) for i=1:3]
        ds = TTest.DRStore(dd)

        @test maximum(cat(1, ds.data...) - ds.flat) == 0.0


        ds.flat[end] = 10
        @test maximum(cat(1, ds.data...) - ds.flat) == 0.0

        ds.data[2][:] *= 0.1
        @test maximum(cat(1, ds.data...) - ds.flat) == 0.0

        @test TTest.distance(ds,ds) == 0.0
        # @time TTest.distance(ds,ds0)

    end

    begin "sarrays"

        ddr1 = DRStore(SVector{2,Float64}, [20000,30000,40000])
        ddr2 = DRStore(SVector{2,Float64}, [20000,30000,40000])

        @time maxabs(ddr1)
        @time distance(ddr1, ddr2)

        # tbc...

    end
end
