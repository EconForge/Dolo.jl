var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Algorithms",
    "title": "Algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Algorithms-1",
    "page": "Algorithms",
    "title": "Algorithms",
    "category": "section",
    "text": "CurrentModule = Dolo"
},

{
    "location": "index.html#Dolo.time_iteration",
    "page": "Algorithms",
    "title": "Dolo.time_iteration",
    "category": "Function",
    "text": "Computes a global solution for a model via backward time iteration. The time iteration is applied to the residuals of the arbitrage equations.\n\nIf the initial guess for the decision rule is not explicitly provided, the initial guess is provided by ConstantDecisionRule. If the stochastic process for the model is not explicitly provided, the process is taken from the default provided by the model object, model.exogenous\n\nArguments\n\nmodel::NumericModel: Model object that describes the current model environment.\nprocess: The stochastic process associated with the exogenous variables in the model.\ninit_dr: Initial guess for the decision rule.\n\nReturns\n\ndr: Solved decision rule.\n\n\n\n"
},

{
    "location": "index.html#Time-iteration-algorithm:-1",
    "page": "Algorithms",
    "title": "Time iteration algorithm:",
    "category": "section",
    "text": "time_iterationGlobalverbose::Bool=true\ndetails::Bool=true\ninterpolation (multilinear, cubic, chebychev, ...)\n...\ndiscretization:\ntype::Symbol=rouwenhorst\noptions: (rouwenhorst)\nN::Int=3\noptions: (tauchen)\nN::Int=3\nmu::Int=3, number of standard deviations\noptions: (gdp)\nmu::Float=3, number of standard deviations\nN::5, number of points today\nN_int::5 number of integration nodes in each dimension\nouter loop control:\nmaxit::Int=100\ntol_ϵ::Float=1e-8\ntol_η::Float=(1-λ)*tol_ϵ\nhow to estimate lambda (def: 0.99)\ncomplementarities:\ncomplementarities::Bool=true\ninner loop control:\nsolver-type= safeguarded-newton\nsolver-options: (for newton)\nverbose::Bool=false\nmaxit::Int=10\ntol::Float=1e-6\neps::Float=1e-8 (to evaluate numerical jacobian)\nn_bsteps::Int=5\nlam_bsteps::Float=0.5\nwhether to use expectation function\nuse_expectations::Bool=false (true requires functions fh)\nuse_direct_response::Bool=false (uses explicit formula to to solve arbitage equations, requires dh)Full signature (for now):time_iteration(model, discretized_process, endogenous_grid, initial_guess;\n    maxit, tol_ϵ, tol_η, λ,\n    complementarities,\n    grid=Dict(...), # overrides get_grid() values\n    solver=Dict(...),  # options for solver, unpacked an passed to solver\n    domain=Dict(...), # overrides get_domain() values\n    use_expectations, use_direct_response\n)"
},

{
    "location": "index.html#Dolo.time_iteration_direct",
    "page": "Algorithms",
    "title": "Dolo.time_iteration_direct",
    "category": "Function",
    "text": "Computes a global solution for a model via backward time iteration. The time iteration is  applied directly to the decision rule of the model.\n\nIf the initial guess for the decision rule is not explicitly provided, the initial guess is provided by ConstantDecisionRule. If the stochastic process for the model is not explicitly provided, the process is taken from the default provided by the model object, model.exogenous.\n\nArguments\n\nmodel::NumericModel: Model object that describes the current model environment.\nprocess: The stochastic process associated with the exogenous variables in the model.\ninit_dr: Initial guess for the decision rule.\n\nReturns\n\ndr: Solved decision rule.\n\n\n\n"
},

{
    "location": "index.html#Time-iteration-algorithm-(direct):-1",
    "page": "Algorithms",
    "title": "Time iteration algorithm (direct):",
    "category": "section",
    "text": "time_iteration_direct"
},

{
    "location": "index.html#Dolo.value_iteration",
    "page": "Algorithms",
    "title": "Dolo.value_iteration",
    "category": "Function",
    "text": "Solve for the value function and associated decision rule using value function iteration.\n\nArguments\n\nmodel::NumericModel: Model object that describes the current model environment.\npdr: Initial guess for the decision rule.\n\nReturns\n\ndr: Solved decision rule object.\ndrv: Solved value function object.\n\n\n\n"
},

{
    "location": "index.html#Value-iteration-1",
    "page": "Algorithms",
    "title": "Value iteration",
    "category": "section",
    "text": "value_iterationOptions:verbose::Bool=true\ndetails::Bool=true\ninterpolation (multilinear, cubic, chebychev, ...)\n...\ndiscretization:\nmethod::Symbol=rouwenhorst\nN::Int=3\nouter loop control:\nmaxit::Int=100\ntol_ϵ::Float=1e-8\ntol_η::Float=(1-λ)*tol_ϵ\nhow to estimate lambda (def: 0.99)\ncomplementarities:\ncomplementarities::Bool=true\ninner loop control:\nsolver-type= safeguarded-newton\nsolver-options:\nverbose::Bool=false\nmaxit::Int=10\ntol::Float=1e-6\neps::Float=1e-8 (to evaluate numerical jacobian)\nn_bsteps::Int=5\nlam_bsteps::Float=0.5\nwhether to use expectation function\nuse_expectations::Bool=false (true requires functions fh)\nuse_direct_response::Bool=false (uses explicit formula to to solve arbitage equations, requires dh)"
},

]}
