const path = Pkg.dir("Dolo", "examples", "models", "dynare", "EA_QUEST3.mod")
const sm = Dolo.load_modfile(path)
const ff = Dolang.FunctionFactory(sm, :dynare)

# functions we will benchmark
regex_create_symbolic_model(path) = Dolo.load_modfile(path)
parse_normalize(sm) = Dolang.normalize(sm.equations[:dynare])
build_function_factory(sm) = Dolang.FunctionFactory(sm, :dynare)
build_levels_func_body(ff) = Dolang.func_body(ff, Dolang.Der{0})
compute_jacobian_derivatives(ff) = Dolang._jacobian_expr_mat(ff)
build_jacobian_func_body(ff) = Dolang.func_body(ff, Dolang.Der{1})
compute_hessian_derivatives(ff) = Dolang._hessian_exprs(ff)
build_hessian_func_body(ff) = Dolang.func_body(ff, Dolang.Der{2})
solve_triangular_system(sm) = Dolo.solve_triangular_system(sm.calibration)

# build benchmark group
suite["eaquest"] = BenchmarkGroup(["eaquest"])

for (f, arg) in [(regex_create_symbolic_model, path),
                 (parse_normalize, sm),
                 (build_function_factory, sm),
                 (build_levels_func_body, ff),
                 (compute_jacobian_derivatives, ff),
                 (build_jacobian_func_body, ff),
                 (compute_hessian_derivatives, ff),
                 (build_hessian_func_body, ff),
                 (solve_triangular_system, sm),
                 ]
    # remove `DoloBenchmarks` module prefix before each function name
    k = symbol(split(string(f), ".")[2])
    suite["eaquest"][k]= @benchmarkable $(f)($(arg))
end
