module DoloTests

using Test
using Dolo, DataStructures

tests = length(ARGS) > 0 ? ARGS : [
                                   "model_types",
                                  # #  "model_import",
                                   "test_calibration",
                                   "test_algos",
                                   "test_perfect_foresight",
                                  # #  "test_features",
                                  # #  "test_linter",
                                  # #  "test_minilang",
                                   "test_decision_rules",
                                   "test_discretize_process"
                                 ]

for t in tests
    include("$(t).jl")
end

end  # module
