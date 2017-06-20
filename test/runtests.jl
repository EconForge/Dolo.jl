module DoloTests

using Dolo, DataStructures, Base.Test

tests = length(ARGS) > 0 ? ARGS : ["model_types",
                                   "model_import",
                                   "test_algos",
                                   "test_perfect_foresight"]

for t in tests
    include("$(t).jl")
end

end  # module
