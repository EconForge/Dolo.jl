module DoloTests

using Dolo, DataStructures, Base.Test

tests = length(ARGS) > 0 ? ARGS : ["model_types",
                                   "model_import",
                                   "test_calibration",
                                   "test_algos",
                                #    "test_perfect_foresight",
                                #    "test_features",
                                   "test_minilang",
                                   "test_decision_rules"]

for t in tests
    include("$(t).jl")
end

end  # module
