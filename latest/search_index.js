var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Dolo.jl",
    "title": "Dolo.jl",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Dolo.jl-1",
    "page": "Dolo.jl",
    "title": "Dolo.jl",
    "category": "section",
    "text": "currentModule = Dolo"
},

{
    "location": "index.html#What-is-dolo-?-1",
    "page": "Dolo.jl",
    "title": "What is dolo ?",
    "category": "section",
    "text": "Dolo is a tool to describe and solve economic models. It provides a simple classification scheme to describe many types of models, allows to write the models as simple text files and compiles these files into efficient Julia objects representing them. It also provides many reference solution algorithms to find the solution of these models under rational expectations.Dolo understand several types of nonlinear models with occasionnally binding constraints (with or without exogenous discrete shocks), as well as local pertubations models, like Dynare. It is a very adequate tool to study zero-lower bound issues, or sudden-stop problems, for instance.Sophisticated solution routines are available: local perturbations, perfect foresight solution, policy iteration, value iteration. Most of these solutions are either parallelized or vectorized. They are written in pure Julia, and can easily be inspected or adapted."
},

{
    "location": "index.html#Installation-1",
    "page": "Dolo.jl",
    "title": "Installation",
    "category": "section",
    "text": "To install latest stable release: Pkg.add(\"Dolo\")."
},

{
    "location": "modeling_language.html#",
    "page": "The dolo language",
    "title": "The dolo language",
    "category": "page",
    "text": ""
},

{
    "location": "modeling_language.html#The-dolo-language-1",
    "page": "The dolo language",
    "title": "The dolo language",
    "category": "section",
    "text": "The easiest way to code a model in dolo consists in using specialized Yaml files also referred to as dolo model files."
},

{
    "location": "modeling_language.html#YAML-format-1",
    "page": "The dolo language",
    "title": "YAML format",
    "category": "section",
    "text": "YAML stands for Yet Another Markup Language. It is a serialization language that allows complex data structures in a human-readable way. Atomic elements are floats, integers and strings. An ordered list can be defined by separating elements with commas and enclosing them with square brackets:[1,2,3]Equivalently, it can be done on several lines, by prepending - to each line- 'element'\n- element         # quotes are optional there is no ambiguity\n- third element   # this is interpreted as ``'third element'``Associative arrays map keys(simple strings) to arbitrary values as in the following example:{age: 18, name: peter}Mappings can also be defined on several lines, and structures can be nested by using indentation (use spaces no tabs):age: 18\nname: peter\noccupations:\n  - school\n  - guitar\nfriends:\n  paula: {age: 18}The correspondance between the yaml definition and the resulting Julia object is very transparent. YAML mappings and lists are converted to Julia dictionaries and arrays respectively.Special objects from the Dolo language can be created by adding a tag to yaml nodes as in the following examples:- !AR1:\n    rho: 0.9\n    Sigma: [[0.1]]or- !Product:\n     !AR1:\n        rho: 0.9\n        Sigma: [[0.1]]\n     !AgingProcess:\n          mu: 0.01     # death probability\n          K: 10        # number of agesAny model file must be syntactically correct in the Yaml sense, before the content is analysed further. More information about the YAML syntax can be found on the YAML website, especially from the language specification."
},

{
    "location": "modeling_language.html#Example-1",
    "page": "The dolo language",
    "title": "Example",
    "category": "section",
    "text": "Here is an example model contained in the file examples\\models\\rbc.yamlThis model can be loaded using the command:model = yaml_import(`examples\\global_models\\example.yaml`)The function yaml_import (cross) will raise errors until the model satisfies basic compliance tests. [more of it below]. In the following subsections, we describe the various syntaxic rules prevailing while writing yaml files."
},

{
    "location": "modeling_language.html#Sections-1",
    "page": "The dolo language",
    "title": "Sections",
    "category": "section",
    "text": "A dolo model consists of the following 4 or 5 parts:a symbols section where all symbols used in the model must be   defined\nan equations section containing the list of equations\na calibration section providing numeric values for the symbols\nan options section containing additional information\na covariances or markov_chain section where exogenous shocks are   definedThese section have context dependent rules. We now review each of them in detail:"
},

{
    "location": "modeling_language.html#Declaration-section-1",
    "page": "The dolo language",
    "title": "Declaration section",
    "category": "section",
    "text": "This section is introduced by the symbols keyword. All symbols appearing in the model must be defined there.Symbols must be valid Julia identifiers (alphanumeric not beginning with a number) and are case sensitive. Greek letters (save for lambda which is a keyword) are recognized. Subscripts and superscripts can be denoted by _ and __ respectively. For instance beta_i_1__d will be printed nicely as beta_i1^d.Symbols are sorted by type as in the following example:symbols:\n  variables: [a, b]\n  shocks: [e]\n  parameters: [rho]Note that each type of symbol is associated with a symbol list (as [a,b]).noteA common mistake consists in forgetting the commas, and using spaces only. This doesn't work since two symbols are recognized as one.The expected types depend on the model that is being written:For Dynare models, all endogenous variables must be listed as   variables with the exogenous shocks being listed as shocks (as in   the example above).noteThe variables, shocks and parameters keywords correspond to the var, varexo and param keywords in Dynare respectively.Global models require the definition of the parameters, and to providea list of states and controls. Mixed states model also require markov_states that follow a discrete markov chain, while continuous states model need to identify the i.i.d shocks that hit the model. If the corresponding equations are given (see next subsection) optional symbols can also be defined. Among them: values, expectations."
},

{
    "location": "modeling_language.html#Declaration-of-equations-1",
    "page": "The dolo language",
    "title": "Declaration of equations",
    "category": "section",
    "text": "The equations section contains blocks of equations sorted by type.Epxressions follow (roughly) the Dynare conventions. Common arithmetic operators (+,-,*,/,\\^) are allowed with conventional priorities as well as usual functions (sqrt, log, exp, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh). The definitions of these functions match the definitions from the numpy package. All symbols appearing in an expression must either be declared in the symbols section or be one of the predefined functions. Any symbol s that is not a parameter is assumed to be considered at date t. Values at date t+1 and t-1 are denoted by s(1) and s(-1) respectively.All equations are implicitly enclosed by the expectation operator E_tleftcdots right. Consequently, the law of motion for the capitalk_t+1 = (1-delta) k_t +  i_t + epsilon_tis written as:k = (1-delta)*k(-1) + i(-1)while the Euler equationE_t left 1=beta left( fracc_t+1c_t + (1-delta)+r_t+1 right) rightis translated by:1 = beta*(c/c(1))^(sigma)*(1-delta+rk(1))An equation can consist of one expression, or two expressions separated by =. There are two types of equation blocks:condition blocks\nIn these blocks, each equation lhs = rhs define the scalar value (rhs)-(lhs). A list of of such equations, i.e a block, defines a multivariate function of the appearing symbols. Certain condition blocks, can be associated with complementarity conditions separated by |` as in rhs-lhs | 0 < x < 1. In this case it is advised to omit the equal sign in order to make it easier to interpret the complementarity. Also, when complementarity conditions are used, the ordering of variables appearing in the complementarities must match the declaration order (more in section Y).\ndefinition blocksDefinition blocks differ from condition blocks in that they define a group of variables (states or auxiliaries) as a function of the right hand side.The types of variables appearing on the right hand side depend on the block type. The variables enumerated on the left hand-side must appear in the declaration order.noteIn the RBC example, the auxiliary block defines variables (y,c,rk,w) that can be directly deduced from the states and the controls:auxiliary:\n    - y = z*k^alpha*n^(1-alpha)\n    - c = y - i\n    - rk = alpha*y/k\n    - w = (1-alpha)*y/wNote that the declaration order matches the order in which variables appear on the left hand side. Also, these variables are defined recursively: c, rk and w depend on the value for y. In contrast to the calibration block, the definition order matters. Assuming that variables were listed as (c,y,rk,w) the following block would provide incorrect result since y is not known when c is evaluated.auxiliary:\n    - c = y - i\n    - y = z*k^alpha*n^(1-alpha)\n    - rk = alpha*y/k\n    - w = (1-alpha)*y/w"
},

{
    "location": "modeling_language.html#Calibration-section-1",
    "page": "The dolo language",
    "title": "Calibration section",
    "category": "section",
    "text": "The role of the calibration section consists in providing values for the parameters and the variables. The calibration of all parameters appearing in the equation is of course strictly necessary while the  calibration of other types of variables is useful to define the steady-state or an initial guess of the steady-state.The calibrated values are also substituted in other sections, including the shocks and options section. This is particularly useful to make the covariance matrix depend on model parameters, or to adapt the state-space to the model's calibration.The calibration is given by an associative dictionary mapping symbols to define with values. The values can be either a scalar or an expression. All symbols are treated in the same way, and values can depend upon each other as long as there is a way to resolve them recursively.In particular, it is possible to define a parameter in order to target a special value of an endogenous variable at the steady-state. This is done in the RBC example where steady-state labour is targeted with n: 0.33 and the parameter phi calibrated so that the optimal labour supply equation holds at the steady-state (chi: w/c^sigma/n^eta).All symbols that are defined in the symbols section but do not appear in the calibration section are initialized with the value nan without issuing any warning.noteNo clear policy has been established yet about how to deal with undeclared symbols in the calibration section. Avoid them."
},

{
    "location": "modeling_language.html#Shock-specification-1",
    "page": "The dolo language",
    "title": "Shock specification",
    "category": "section",
    "text": "The way shocks are specified depends on the type of model. They are constructed using a the rules for mini-languages defined in section [ref]."
},

{
    "location": "modeling_language.html#Distribution-1",
    "page": "The dolo language",
    "title": "Distribution",
    "category": "section",
    "text": "For Dynare and continuous-states models, one has to specifiy a multivariate distribution of the i.i.d. process for the vector of shocks (otherwise shocks are assumed to be constantly 0). This is done in the distribution section. A gaussian distrubution (only one supported so far), is specified by supplying the covariance matrix as a list of lists as in the following example.distribution:\n\n    Normal: [\n            [sigma_1, 0.0],\n            [0.0, sigma_2]\n        ]"
},

{
    "location": "modeling_language.html#Markov-chains-1",
    "page": "The dolo language",
    "title": "Markov chains",
    "category": "section",
    "text": "When the model is driven by an exogenous discrete markov chain, that is for DTMSCC models, shocks are defined in the discrete_transition section. The objects allowed in this section are: MarkovChain, AR1, MarkovTensormarkov chain can be constructed in several ways:by listing directly a list of states, and a transition matrix as   in :\ndiscrete_transition:\n    MarkovChain:   # a markov chain is defined by providing:\n        - [ [0.0, -0.02]           # a list of markov states\n            [0.0,  0.02]\n            [-0.1, 0.02]]\n        - [ [ 0.98, 0.01, 0.01],   # a transition matrix\n            [ 0.10, 0.01, 0.90],\n            [ 0.05, 0.05, 0.90] ]by using primitives to construct a discretized process from an   AR1:\ndiscrete_transition:\n    AR1:\n        rho: 0.9\n        sigma: [\n                [0.01, 0.001]\n                [0.001, 0.02]\n            ]\n        N: 3\n        method: rouwenhorst   # the alternative is tauchen\nby combining two processes together:\ndiscrete_transition:\n    MarkovTensor:\n        - AR1:\n            rho: 0.9\n            sigma: [\n                    [0.01, 0.001]\n                    [0.001, 0.02]\n                ]\n            N: 3\n            method: rouwenhorst   # the alternative is tauchen\n        - AR1:\n            rho: 0.9\n            sigma: 0.01\n            N: 2\n            method: rouwenhorst   # the alternative is tauchen"
},

{
    "location": "modeling_language.html#Domain-1",
    "page": "The dolo language",
    "title": "Domain",
    "category": "section",
    "text": "The domain section defines lower and upper bounds for the exogenous and endogenous states. For example, in the RBC model, we write:domain:\n  z: [-2*sig_z/(1-rho^2)^0.5,  2*sig_z/(1-rho^2)^0.5]\n  k: [ k*0.5, k*1.5]The part for z sets the bounds for the productivity process to be two times its asymptotic standard deviation. The boundaries for capital are a 50% bracket around its steady-state level. "
},

{
    "location": "modeling_language.html#Options-1",
    "page": "The dolo language",
    "title": "Options",
    "category": "section",
    "text": "The options sections contains extra information needed to solve the model.The section follows the mini-language convention, with all calibrated values replaced by scalars and all keywords allowed.Here we can define the grids and specify their type, Cartesian. We can also specify how many grid points we want for z and k. Here we choose 5 points for z and 50 points for k. Note that the grid points are listed in accordance with the declaration order of the variables.options:\n  grid: !Cartesian\n        orders: [5, 50]"
},

{
    "location": "model_specification.html#",
    "page": "Model Specification",
    "title": "Model Specification",
    "category": "page",
    "text": ""
},

{
    "location": "model_specification.html#Model-Specification-1",
    "page": "Model Specification",
    "title": "Model Specification",
    "category": "section",
    "text": ""
},

{
    "location": "model_specification.html#Variable-types-1",
    "page": "Model Specification",
    "title": "Variable types",
    "category": "section",
    "text": "The following types of variables can be used in models:exogenous (m)\nstates (s)\ncontrols (x)\nauxiliaries (y)\nrewards (r)\nvalues (v)\nexpectations (z)\nparameters (p)Symbol types that are present in a model are always listed in that order."
},

{
    "location": "model_specification.html#State-space-1",
    "page": "Model Specification",
    "title": "State-space",
    "category": "section",
    "text": "Decisions are characterized by a vector m of exogenous variables (exogenous states) and by a vector s of endogenous states. The unknown vector of controls x is a function varphi of the states such that:x = varphi(ms)The function varphi is typically approximated by the solution algorithm. It can be either a Taylor expansion, or an interpolating object (splines, smolyak). Once obtained, it can be evaluated at any point of the the state space (represented by a couple of vectors), or at a list of points (of couple of matrices):dr = time_iteration(model)\nm0 = model.calibration[:states]\ns0 = model.calibration[:states]\ndr(m0,s0)                               # evaluates at the steady-state\ndr([0.0;-0.01;-0.1], [2.5;2.5;2.5])     # evaluates at a list of pointsA decision rule is defined on two discretized grids: one for the exogenous states and one for the endogenous ones. If the exogenous grid points are numbered by i, we can define non-ambiguously x = varphi(is) as  varphi(m_is), which corresponds to the following code:dr(2,s0)                 # evaluates at the second exogenous point\ndr(2, [2.5;2.5;2.5])     # evaluates at a list of points"
},

{
    "location": "model_specification.html#Valid-equations-1",
    "page": "Model Specification",
    "title": "Valid equations",
    "category": "section",
    "text": "The various equations understood by Dolo are descrbed below:<!– would be cool to know how to make comments ;-) –>"
},

{
    "location": "model_specification.html#Transitions-1",
    "page": "Model Specification",
    "title": "Transitions",
    "category": "section",
    "text": "- name: `transition`\n- short name: `g`Transitions are given by a function g such that at all times:s_t = g(m_t-1 s_t-1 x_t-1 m_t)where m_t is a vector-valued exogenous process.noteIn the RBC model, the vector of endogenous states is s_t=(a_tk_t). The transitions are:a_t = rho a_t-1 + epsilon_tk_t = (1-delta)k_t-1 + i_t-1If epsilon_t follow a normal distribution with variance sigma_epsilon, the yaml file is amended with:symbols:\n    states: [a,k]\n    controls: [i]\n    exogenous: [epsilon]\n    ...\nequations:\n    transition:\n        a = rho*a(-1) + e\n        k = k(-1)*(1-delta) + i(-1)\nexogenous: !Normal\n    Sigma: [[sigma_epsilon]]Note that transition equations must list states in declaration order. Also, in this example the productivity process a_t is essentially an exogenous AR1 process. Declaring it as such (instead of epsilon), provides additional information to the solvers and can lead to faster solution time."
},

{
    "location": "model_specification.html#Auxiliary-variables-/-Definitions-1",
    "page": "Model Specification",
    "title": "Auxiliary variables / Definitions",
    "category": "section",
    "text": "- name: `auxiliary`\n- short name: `a`In order to reduce the number of variables, it is useful to define auxiliary variables y_t using a function a such that:y_t = a(m_t s_t x_t)These variables are defined in a special definitions block, outside of equations. When auxiliary variables appear in an equation they are automatically substituted by the corresponding expression in m_t,s_t and x_t.noteIn the RBC model, three auxiliary variables are defined y_t c_t r_kt and w_t. They are a closed form function of a_t k_t i_t n_t. Defining these variables speeds up computation since they are don't need to be solved for or interpolated."
},

{
    "location": "model_specification.html#Utility-function-and-Bellman-equation-1",
    "page": "Model Specification",
    "title": "Utility function and Bellman equation",
    "category": "section",
    "text": "- name: `utility`\n- short name: `u`The (separable) value equation defines the value v_t of a given policy as:v_t = u(m_ts_tx_t) + beta E_t left v_t+1 rightThis gives rise to the Bellman eqution:v_t = max_x_t left( u(m_t s_tx_t) + beta E_t left v_t+1 right right)These two equations are characterized by the reward function u and the discount rate beta. Function u defines the vector of symbols rewards. Since the definition of u alone is not sufficient, the parameter used for the discount factor must be given to routines that compute the value. Several values can be computed at once, if U is a vector function and beta a vector of discount factors, but in that case in cannot be used to solve for the Bellman equation.noteOur RBC example defines the value as v_t = frac(c_t)^1-gamma1-gamma-chi frac(n_t)^1+eta1+eta + beta E_t v_t+1. This information is coded using: symbols:\n    ...\n    rewards: [r]\n\nequations:\n    ...\n    utility:\n        - r = c^(1-gamma)/(1-gamma)- chi*n^(1+eta)/(1+eta)\n\ncalibration:\n    ...\n    beta: 0.96   # beta is the default name of the discount"
},

{
    "location": "model_specification.html#Value-1",
    "page": "Model Specification",
    "title": "Value",
    "category": "section",
    "text": "- name: `value`\n- short name: `w`A more general updating equation can be useful to express non-separable utilities or prices. The vector of (generalized) values v^* are defined by a function w such that:v_t = w(m_ts_tx_tv_tm_t+1s_t+1x_t+1v_t+1)As in the separable case, this function can either be used to compute the value of a given policy x=varphi() or in order solve the generalized Bellman equation:v_t = max_x_t left( w(m_ts_tx_tv_tm_t+1s_t+1x_t+1v_t+1) right)noteInstead of defining the rewards of the RBC example, one can instead define a value updating equation:symbols:\n    ...\n    values: [v]\n\nequations:\n    ...\n    value:\n        - v = c^(1-gamma)/(1-gamma)*(1-n...) + beta*v(1)"
},

{
    "location": "model_specification.html#Boundaries-1",
    "page": "Model Specification",
    "title": "Boundaries",
    "category": "section",
    "text": "- name: `controls_lb` and `controls_ub`\n- short name: `lb` and `ub`The optimal controls must also satisfy bounds that are function of states. There are two functions underlineb() and overlineb() such that:underlineb(m_t s_t) leq x_t leq overlineb(m_t s_t)noteIn our formulation of the RBC model we have excluded negative investment, implying i_t geq 0. On the other hand, labour cannot be negative so that we add lower bounds to the model:equations:\n    ...\n    controls_lb:\n        i = 0\n        n = 0Specifying the lower bound on labour actually has no effect since agents endogeneously choose to work a positive amount of time in order to produce some consumption goods. As for upper bounds, it is not necessary to impose some: the maximum amount of investment is limited by the Inada conditions on consumption. As for labour n, it can be arbitrarly large without creating any paradox. Thus the upper bounds are omitted (and internally treated as infinite values)."
},

{
    "location": "model_specification.html#Euler-equation-1",
    "page": "Model Specification",
    "title": "Euler equation",
    "category": "section",
    "text": "- name: `arbitrage` (`equilibrium`)\n- short name: `f`A general formulation of the Euler equation is:0 = E_t left f(m_t s_t x_t m_t+1 s_t+1 x_t+1) rightNote that the Euler equation and the boundaries interact via \"complentarity equations\". Evaluated at one given state, with the vector of controls x=(x_1  x_i  x_n_x), the Euler equation gives us the residuals r=(f_1  f_i  f_n_x). Suppose that the i-th control x_i is supposed to lie in the interval  underlineb_i overlineb_i . Then one of the following conditions must be true:f_i = 0\nf_i0 and x_i=overlineb_i\nf_i0 and x_i=underlineb_iBy definition, this set of conditions is denoted by:f_i = 0 perp underlineb_i leq x_i leq overlineb_iThese notations extend to a vector setting so that the Euler equations can also be written:0 = E_t left f(m_t s_t x_t m_t+1 s_t+1 x_t+1) right perp underlineb(m_t s_t) leq x_t leq overlineb(m_t s_t)Specifying the boundaries together with Euler equation, or providing them separately is exactly equivalent. In any case, when the boundaries are finite and occasionally binding, some attention should be devoted to write the Euler equations in a consistent manner. In particular, note that the Euler equations are order-sensitive.The Euler conditions, together with the complementarity conditions typically often come from Kuhn-Tucker conditions associated with the Bellman problem, but that is not true in general.noteThe RBC model has two Euler equations associated with investment and labour supply respectively. They are added to the model as:arbitrage:\n    - 1 - beta*(c/c(1))^(sigma)*(1-delta+rk(1))   | 0 <= i <= inf\n    - w - chi*n^eta*c^sigma                       | 0 <= n <= infPutting the complementarity conditions next to the Euler equations, instead of entering them as separate equations, helps to check the sign of the Euler residuals when constraints are binding. Here, when investment is less desirable, the first expression becomes bigger. When the representative agent is prevented from investing less due to the constraint (i.e. i_t=0), the expression is then positive, consistent with the complementarity conventions."
},

{
    "location": "model_specification.html#Expectations-1",
    "page": "Model Specification",
    "title": "Expectations",
    "category": "section",
    "text": "- name: `expectation`\n- short name: `h`The vector of explicit expectations z_t is defined by a function h such that:z_t = E_t left h(m_t+1 s_t+1x_t+1) rightIn the RBC example, one can define the expected value tomorrow of one additional unit invested tomorrow:\n\n.. math::\n\n    m_t=\\beta*(c_{t+1}^(-\\sigma)*(1-\\delta+r_{k,t+1})\n\n It is a pure expectational variable in the sense that it is solely determined by future states and decisions. In the model file, it would be defined as:\n\n.. code: yaml\n\n    symbols:\n        ...\n        expectations: [z]\n\n    equations:\n        ...\n        - z = beta*(c(1))^(-sigma)*(1-delta+rk(1))"
},

{
    "location": "model_specification.html#Generalized-expectations-1",
    "page": "Model Specification",
    "title": "Generalized expectations",
    "category": "section",
    "text": "- name: `expectation_2`\n- short name: `h_2`The vector of generalized explicit expectations z_t is defined by a function h^star such that:z_t = E_t left h^star(m_ts_tx_tm_t+1s_t+1x_t+1) right"
},

{
    "location": "model_specification.html#Euler-equation-with-expectations-1",
    "page": "Model Specification",
    "title": "Euler equation with expectations",
    "category": "section",
    "text": "- name: `arbitrage_2` (`equilibrium_2`)\n- short name: `f_2`If expectations are defined using one of the two preceding definitions, the Euler equation can be rewritten as:0 = f(m_t s_t x_t z_t) perp underlineb(m_t s_t) leq x_t leq overlineb(m_t s_t)noteGiven the definition of the expectation variable m_t, today's consumption is given by: c_t = z_t^left(-frac1sigmaright) so the Euler equations are rewritten as:arbitrage_2:\n    - 1 - beta*(c)^(sigma)/m   | 0 <= i <= inf\n    - w - chi*n^eta*c^sigma    | 0 <= n <= infNote the type of the arbitrage equation (arbitrage_2 instead of arbitrage).However c_t is not a control itself,but the controls i_t n_t can be easily deduced:..math:n_t = ((1-\\alpha) z_t k_t^\\alpha m_t/chi)^(1/(eta+\\alpha))\ni_t = z_t k_t^\\alpha n_t^(1-\\alpha) - (m_t)^(-1/\\sigma)This translates into the following YAML code:equations:\n    - n = ((1-alpha)*a*k^alpha*m/chi)^(1/(eta+alpha))\n    - i = z*k^alpha*n^(1-alpha) - m^(-1/sigma)"
},

{
    "location": "model_specification.html#Direct-response-function-1",
    "page": "Model Specification",
    "title": "Direct response function",
    "category": "section",
    "text": "- name: `direct_response`\n- short name: `d`In some simple cases, there a function d() giving an explicit definition of the controls:x_t = d(m_t s_t z_t)Compared to the preceding Euler equation, this formulation saves computational time by removing the need to solve a nonlinear system to recover the controls implicitly defined by the Euler equation."
},

{
    "location": "algos.html#",
    "page": "Solution Algorithms",
    "title": "Solution Algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "algos.html#Solution-Algorithms-1",
    "page": "Solution Algorithms",
    "title": "Solution Algorithms",
    "category": "section",
    "text": ""
},

{
    "location": "algos.html#Dolo.value_iteration",
    "page": "Solution Algorithms",
    "title": "Dolo.value_iteration",
    "category": "Function",
    "text": "Solve for the value function and associated decision rule using value function iteration.\n\nArguments\n\nmodel::NumericModel: Model object that describes the current model environment.\npdr: Initial guess for the decision rule.\n\nReturns\n\ndr: Solved decision rule object.\ndrv: Solved value function object.\n\n\n\n"
},

{
    "location": "algos.html#Value-iteration-1",
    "page": "Solution Algorithms",
    "title": "Value iteration",
    "category": "section",
    "text": "Value function algorithms require a Bellman representation of the model:s_t = colorred g left( m_t-1 s_t-1 x_t-1 m_t right)V(m_t s_t) = max_underlinex(m_ts_t) leq x_t leq overlinex(m_ts_t) colorred uleft(m_t s_t x_tright) + colorblue beta E_m_t Vleft(m_t+1 s_t+1right)where g is the transition function, u the instantaneous felicity function and beta the time-discount parameter.The solution of this problems produces naturally to functions x=varphi(ms) and  V=varphi(ms) for the controls and the value function respectively.Given an initial value function V^n(ms) applying the maximum operator produces a new, improved, value function tildeV^n(ms) and a corresponding policy rule x=varphi^n(ms). This is an improvement step.Note that at this stage, tildeV^n is not the value of following varphi^n forever. An evaluation step performs the recursion (note the absence of a max):x_t=varphi^n(m_s s_t)s_t = colorred g left( m_t-1 s_t-1 x_t-1 m_t right)V^k+1(m_t s_t) = colorred uleft(m_t s_t x_tright) + colorblue beta E_m_t V^kleft(m_t+1 s_t+1right)Starting from the initial guess, tildeV^n+1, this recursion converges to the value V^n+1() associated to varphi^n+1(). If high accuracy is not required, it is common to restrict the number of steps to perform a partial evaluation.The value function algorithm implemented in Dolo follows the following scheme:Given an initial guess for the policy rule varphi^0(), evaluate the corresponding value function V^0().\nThen given varphi_n:\nperform one improvement step to get varphi^n+1(), tildeV^n+1()\nperform a partial evaluation of varphi^n+1()\nstarting from tildeV^n+1() to get V^n+1() using at most K steps.\nCompute eta_n+1=varphi^n-varphi^n+1 and epsilon_n+1=varphi^n-varphi^n+1\nif eta_n+1tau_eta and epsilon_n+1tau_epsilon, return\nelse return to step 2.This algorihtm is controlled by the precision parameters tau_eta, tau_epsilon and the evaluation length K. Note that setting K=0 corresponds to the naive (but robust ?) VFI algorithm while K high corresponds to Howard improvements which converge faster: convergence of the outer loop varphi^n is quadratic instead of geometric.value_iteration"
},

{
    "location": "algos.html#Dolo.time_iteration",
    "page": "Solution Algorithms",
    "title": "Dolo.time_iteration",
    "category": "Function",
    "text": "Computes a global solution for a model via backward time iteration. The time iteration is applied to the residuals of the arbitrage equations.\n\nIf the initial guess for the decision rule is not explicitly provided, the initial guess is provided by ConstantDecisionRule. If the stochastic process for the model is not explicitly provided, the process is taken from the default provided by the model object, model.exogenous\n\nArguments\n\nmodel::NumericModel: Model object that describes the current model environment.\nprocess: The stochastic process associated with the exogenous variables in the model.\ninit_dr: Initial guess for the decision rule.\n\nReturns\n\ndr: Solved decision rule.\n\n\n\n"
},

{
    "location": "algos.html#Dolo.time_iteration_direct",
    "page": "Solution Algorithms",
    "title": "Dolo.time_iteration_direct",
    "category": "Function",
    "text": "Computes a global solution for a model via backward time iteration. The time iteration is  applied directly to the decision rule of the model.\n\nIf the initial guess for the decision rule is not explicitly provided, the initial guess is provided by ConstantDecisionRule. If the stochastic process for the model is not explicitly provided, the process is taken from the default provided by the model object, model.exogenous.\n\nArguments\n\nmodel::NumericModel: Model object that describes the current model environment.\nprocess: The stochastic process associated with the exogenous variables in the model.\ninit_dr: Initial guess for the decision rule.\n\nReturns\n\ndr: Solved decision rule.\n\n\n\n"
},

{
    "location": "algos.html#Time-iteration-1",
    "page": "Solution Algorithms",
    "title": "Time iteration",
    "category": "section",
    "text": "We consider a model with the form:s_t = gleft(m_t-1 s_t-1 x_t-1 m_t right)0 = E_t left fleft(m_t s_t x_t m_t+1 s_t+1 x_t+1 right) rightwhere g is the state transition function, and f is the arbitrage equation.The time iteration algorithm consists in approximating the optimal controls as a function of exogenous and endogenous controls x_t = varphi(m_ts_t). At step n, the current guess for the control, x(s_t) = varphi^n(m_t s_t), serves as the control being exercised next period.Here is an outline of the algorithm:Start with initial guess varphi^0\nGiven current guess, find the current period's  varphi^n+1(m_ts_t) controls for any (m_ts_t) by solving (numerically)  the arbitrage equation :0 = E_t left fleft(m_t s_t varphi^n+1(m_t s_t) g(s_t varphi^n+1(m_t s_t)) varphi^n(m_t+1g(s_t varphi^n+1(s_t))) right) rightCompute successive approximation errors eta_n=varphi^n+1-varphi^n.\nif eta_n smaller thatn criterion epsilon_eta, return\notherwise return to step 2time_iterationIn some cases, the solution of the Euler equation, can be obtained faster if a direct solution for optimal controls is known as a function expectation as in the following specification:s_t = gleft(m_t-1 s_t-1 x_t-1 m_t right)z_t = E_t left hleft(m_t+1 s_t+1 x_t+1 right) rightx_t = d(m_t s_t z_t)This information can be used by passing the solver=Dict(:type=>:direct) option to the time_iteration function, or by using the devoted function:time_iteration_direct"
},

{
    "location": "algos.html#Dolo.improved_time_iteration",
    "page": "Solution Algorithms",
    "title": "Dolo.improved_time_iteration",
    "category": "Function",
    "text": "Computes a global solution for a model via backward Improved Time Iteration. The algorithm is applied to the residuals of the arbitrage equations. The idea is to solve the system G(x) = 0 as a big nonlinear system in x, where the inverted Jacobian matrix is approximated by an infinite sum (Neumann series).\n\nIf the initial guess for the decision rule is not explicitly provided, the initial guess is provided by ConstantDecisionRule. If the stochastic process for the model is not explicitly provided, the process is taken from the default provided by the model object, model.exogenous\n\nArguments\n\nmodel::NumericModel: Model object that describes the current model environment.\ndprocess: The stochastic process associated with the exogenous variables in the model.\ninit_dr: Initial guess for the decision rule.\nmaxbsteps Maximum number of backsteps.\nverbose Set \"true\" if you would like to see the details of the infinite sum convergence.\nsmaxit Maximum number of iterations to compute the Neumann series.\ncomplementarities\ncompute_radius\ndetails If false returns only a decision rule dr\n\nReturns\n\ndr: Solved decision rule.\ndetails about the iterations is specified.\n\n\n\n"
},

{
    "location": "algos.html#Improved-Time-iteration-1",
    "page": "Solution Algorithms",
    "title": "Improved Time iteration",
    "category": "section",
    "text": "improved_time_iteration"
},

{
    "location": "algos.html#Dolo.perfect_foresight",
    "page": "Solution Algorithms",
    "title": "Dolo.perfect_foresight",
    "category": "Function",
    "text": "Document pf.\n\n\n\n"
},

{
    "location": "algos.html#Perfect-Foresight-1",
    "page": "Solution Algorithms",
    "title": "Perfect Foresight",
    "category": "section",
    "text": "Consider a series for the exogenous process (m_t)_0 leq t leq T. The perfect foresight problem consists in finding the path of optimal controls (x_t)_0 leq t leq T and corresponding states (s_t)_0 leq t leq T such that:s_t = gleft(m_t-1 s_t-1 x_t-1 m_t right)0 = E_t left( fleft(m_t s_t x_t m_t+1 s_t+1 x_t+1right) right)  perp  underlineu = x_t = overlineuSpecial conditions apply for the initial state and controls. Initial state s_0 is given exogenously, or determined so that it corresponds for a steady-state corresponding to m_0. Final states and controls are determined by assuming the exogenous process satisfies m_t=m_T for all tleq T and optimality conditions are satisfied in the last period:f(m_T s_T x_T m_Ts_T x_T) perp underlineu = x_T = overlineu.We assume that underlineu and overlineu are constants. This is not a big restriction since the model can always be reformulated in order to meet that constraint, by adding more equations.The stacked system of equations satisfied by the solution is:Transition Arbitrage\ns_0 = overlines_0 f(m_0 s_0 x_0 m_1 s_1 x_1) perp underlineu = x_0 = overlineu\ns_1 = g(m_0 s_0 x_0 m_1) f(m_1 s_1 x_1 m_2 s_2 x_2) perp underlineu = x_1 = overlineu\n... ...\ns_T = g(m_T-1 s_T-1 x_T-1 m_T) f(m_T s_T x_T m_T s_T x_T) perp underlineu = x_T = overlineuThe system is solved using a nonlinear solver.perfect_foresight"
},

{
    "location": "algos.html#Dolo.residuals",
    "page": "Solution Algorithms",
    "title": "Dolo.residuals",
    "category": "Function",
    "text": "residuals(model::AModel, [calib::ModelCalibration])::Dict\n\nCompute the steady state residuals for the aribtrage and transition equations of model, when these functions are evaluated at the data in calib. If no calib is provided, model.calibration will be used.\n\nSee the docstring for find_deterministic_equilibrium for more information\n\n\n\n"
},

{
    "location": "algos.html#Dolo.find_deterministic_equilibrium",
    "page": "Solution Algorithms",
    "title": "Dolo.find_deterministic_equilibrium",
    "category": "Function",
    "text": "find_deterministic_equilibrium(model::AModel, [calib::ModelCalibration])\n\nSolve for the steady state equilibrium of model, data in cailb to fill in parameter values and provide an initial guess for the states and controls. When no calibration is passed model.calibration is used\n\nThe exogenous variables at time t-1 (m) and t (M) are set to calib[:exogenous].\n\nThe deterministic equilibrium is found by solving for vectors s and x, such that\n\ns = transition(m, s, x, m, p)\n0 = arbitrage(m, s, x, m, s, x, p)\n\n\n\n"
},

{
    "location": "algos.html#Dolo.perturbate",
    "page": "Solution Algorithms",
    "title": "Dolo.perturbate",
    "category": "Function",
    "text": "TBD\n\n\n\n"
},

{
    "location": "algos.html#Local-Analysis-1",
    "page": "Solution Algorithms",
    "title": "Local Analysis",
    "category": "section",
    "text": "residualsfind_deterministic_equilibriumperturbate"
},

{
    "location": "simulate.html#",
    "page": "Inspecting Solutions",
    "title": "Inspecting Solutions",
    "category": "page",
    "text": ""
},

{
    "location": "simulate.html#Inspecting-Solutions-1",
    "page": "Inspecting Solutions",
    "title": "Inspecting Solutions",
    "category": "section",
    "text": ""
},

{
    "location": "simulate.html#Tabulate-Decision-Rule-1",
    "page": "Inspecting Solutions",
    "title": "Tabulate Decision Rule",
    "category": "section",
    "text": "tabulate"
},

{
    "location": "simulate.html#Dolo.response",
    "page": "Inspecting Solutions",
    "title": "Dolo.response",
    "category": "Function",
    "text": "Function \"response\" computes the impulse response functions with several major options:\n\nthe user can provide a vector with the first values of the model's exogenous processes, e1.\nthe user can provide a name of the shock of interest and the size of the shock_name.\nthe user can provide only a name of the shock of interest. The size of the shock is assumed to be a one standard deviation given in the yaml file.\n\nArguments\n\nmodel::NumericModel: Model object that describes the current model environment.\ndr: Solved decision rule.\ne1::ListOfPoints: List of initial model's exogenous processes values.\nIf e1 is not provided, then:\nshock_name: the name of the shock of interest.\noptional:\nImpulse: the size of the shock; default: one standard deviation.\ns0::ListOfPoints: List of initial state variable values; default: model.calibration[:states]\n\nReturns\n\nresponse: Impulse response function.\n\n\n\n"
},

{
    "location": "simulate.html#Impulse-Response-Functions-1",
    "page": "Inspecting Solutions",
    "title": "Impulse Response Functions",
    "category": "section",
    "text": "response"
},

{
    "location": "simulate.html#QuantEcon.simulate",
    "page": "Inspecting Solutions",
    "title": "QuantEcon.simulate",
    "category": "Function",
    "text": "This is the one we document.\n\n\n\n"
},

{
    "location": "simulate.html#Stochastic-Simulation-1",
    "page": "Inspecting Solutions",
    "title": "Stochastic Simulation",
    "category": "section",
    "text": "simulate"
},

]}
