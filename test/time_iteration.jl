root_dir = pkgdir(Dolo)

@testset "Benchmark Models" verbose=true begin

    model_files = ["rbc_iid", "rbc_mc", "rbc_ar1"]

    for fn in model_files

        @testset "$fn" verbose=true begin

            @testset "Julia" verbose=true begin
            
                model = include("$(root_dir)/examples/ymodels/$fn.jl")
                @test isbits(model)
                sol = Dolo.time_iteration(model; verbose=false)

            end


            @testset "YAML" verbose=true begin
            
                model = include("$(root_dir)/examples/ymodels/$fn.jl")
                @test isbits(model)
                sol = Dolo.time_iteration(model; verbose=false)

            end

        end


    end

end
