using Dolang: sanitize
using Lerche: Tree, Token

# used for constructing appropraite dict from YAML object.
function construct_type_map(t::Symbol, constructor::YAML.Constructor,
                            node::YAML.MappingNode)
    mapping = _symbol_dict(YAML.construct_mapping(constructor, node))
    mapping[:tag] = t
    mapping
end

function construct_type_map(t::Symbol, constructor::YAML.Constructor,
                            node::YAML.SequenceNode)
    YAML.construct_sequence(constructor, node)
end

const yaml_types = let
    pairs = [("!Cartesian", :Cartesian),
            ("!Smolyak", :Smolyak),
            ("!Random", :Random),
             ("!Normal", :Normal),
             ("!MarkovChain", :MarkovChain),
             ("!Product", :Product),
             ("!PoissonPro
             function get_symbols(model::Model)
                 syms = model.data[:symbols]
                 symbols = OrderedDict{Symbol,Vector{Symbol}}()
                 for k in keys(syms)
                     symbols[Symbol(k)] = [Symbol(e.value) for e in syms[k]]
                 end
                 return symbols
             endcess", :PoissonProcess),
             ("!DeathProcess", :DeathProcess),
             ("!AgingProcess", :AgingProcess),
             ("!VAR1", :VAR1)]
    Dict{String,Function}([(t, (c, n) -> construct_type_map(s, c, n))
                           for (t, s) in pairs])
end


# mutable struct SModel{ID} <: ASModel{ID}
#     data::Dict{Symbol,Any}
# end

# function SModel(url::AbstractString)
#     if match(r"(http|https):.*", url) != nothing
#         res = HTTP.request("GET", url)
#         txt = (String(res.body))
#         data = _symbol_dict(YAML.load(txt, yaml_types))
#     else
#         data = _symbol_dict(YAML.load_file(url, yaml_types))
#     end
#     id = gensym(basename(url))
#     return SModel{id}(_symbol_dict(data))
# end

# function get_symbols(model::ASModel)
#     syms = model.data[:symbols]
#     symbols = OrderedDict{Symbol,Vector{Symbol}}()
#     for k in keys(syms)
#         symbols[Symbol(k)] = [Symbol(e) for e in syms[k]]
#     end
#     return symbols
# end

# function get_variables(model::ASModel)
#     symbols = get_symbols(model)
#     vars = cat(values(symbols)...; dims=1)
#     dynvars = setdiff(vars, symbols[:parameters] )
# end

# function get_name(model::ASModel)
#     get(model.data, :name, "modeldoesnotwork")
# end

# ### TODO: look for duplicates
# eq_to_expr(eq::Symbol) = eq
# function eq_to_expr(eq::Expr)
#    if !(eq.head==:call && eq.args[1]==:(==))
#       eq
#    else
#       rhs = eq.args[3]
#       lhs = eq.args[2]
#       :($rhs-$lhs)
#    end
# end
# eq_to_expr(eq::Expression) = eq


# function get_equation_block(model, eqname)
#     variables = get_variables(model)
#     variables_str = String[string(e) for e in variables]
#     data = model.rdata
#     yml_node = data["equations"][eqname]

#     if "tag:yaml.org,2002:str" == yml_node.tag
#         block = Dolang.parse_equation_block(yml_node.value; variables=variables_str)
#         equations = block.children
#     else
#         equations = [Dolang.parse_equation(line.value; variables=variables_str) for line in yml_node]
#     end

#     eqs = Union{Expr, Number, Symbol}[]
#     eqs_lb = Union{Expr, Number, Symbol}[]
#     eqs_ub = Union{Expr, Number, Symbol}[]    
#     for eq in equations
#         if eq.data == "double_complementarity"
#             inegality = eq.children[2]
#             # @assert inegality.data == "double_inequality"
#             eq_lb = Dolang.convert(Expr, inegality.children[1])
#             eq_cc = Dolang.convert(Expr, inegality.children[2])
#             eq_ub = Dolang.convert(Expr, inegality.children[3])
#             eq = eq.children[1]
#             push!(eqs_lb, eq_lb)
#             push!(eqs_ub, eq_ub)
#         end
#         if eq.data == "equality"
#             eq = Tree("sub", [eq.children[2], eq.children[1]])          
#         end
#         push!(eqs, Dolang.convert(Expr, eq))
#     end
#     return (eqs, eqs_lb, eqs_ub)
# end



# function get_equations(model::ASModel)

#     eqs = model.data[:equations]
#     recipe = RECIPES[:dtcc]
#     _symbols = get_symbols(model)

#     # prep equations: parse to Expr
#     _eqs = OrderedDict{Symbol,Vector{Expression}}()
#     # for k in Dolo.keys(recipe[:specs])
#     for k in keys(eqs)

#         if !(k in keys(recipe[:specs]))
#             ii=0
#             dist = zeros(length(keys(recipe[:specs])))
#             for kk in keys(recipe[:specs])
#                 ii+=1
#                 dist[ii]=compare(Hamming(), string(kk), string(k))
#             end
#             kk = collect(keys(keys(recipe[:specs]).dict))[findmax(dist)[2]]
#             msg = string("You are using a wrong name for a sub-section title.`$(k)` is not accepted. ",
#                         "Maybe you meant `$(kk)`?")
#             error(msg)
#         end

#         eqs, eq_ub, eq_lb = get_equation_block(model, string(k))

#         _eqs[k] = eqs

#         if k == :arbitrage && (length(eq_ub)>0)
#                 _eqs[:controls_ub] = eq_ub
#                 _eqs[:controls_lb] = eq_lb
#         end

#     end

#     defs = get_definitions(model)

#     dynvars = cat(get_variables(model), [k[1] for k in keys(defs)]...; dims=1)
#     for eqtype in keys(_eqs)
#         _eqs[eqtype] = [sanitize(eq; variables=dynvars) for eq in _eqs[eqtype]]
#     end

#     return _eqs
# end


# function get_infos(model::ASModel)
#     get(model.data, :infos, Dict())
# end

# function get_options(model::ASModel; options=Dict())
#     opts = get(model.data, :options, Dict())
#     return rmerge(opts, options)
# end

# function get_definitions(model::ASModel)::OrderedDict{Tuple{Symbol,Int},SymExpr}
#     # parse defs so values are Expr
#     defs = get(model.data,:definitions, Dict())
#     dynvars = cat(get_variables(model), keys(defs)...; dims=1)
#     _defs = OrderedDict{Tuple{Symbol,Int},SymExpr}(
#         [ (k,0) => Dolang.parse_string(v; variables= dynvars) for (k, v) in defs]
#     )
#     return _defs
# end

# function get_calibration(model::ASModel)
#     calib = get(model.data, :calibration, Dict())
#     # prep calib: parse to Expr, Symbol, or Number
#     _calib  = OrderedDict{Symbol,Union{Expr,Symbol,Number}}()
#     # add calibration for definitions
#     for k in keys(calib)
#         _calib[Symbol(k)] = _expr_or_number(calib[k])
#     end
#     # so far _calib is a symbolic calibration
#     calibration = solve_triangular_system(_calib)
#     symbols = get_symbols(model)
#     return ModelCalibration(calibration, symbols)
# end

# function get_domain(model::ASModel)
#     domain = deepcopy(get(model.data, :domain, Dict()))
#     calib = get_calibration(model)
#     # TODO deal with empty dict and make robust construction
#     states = get_symbols(model)[:states]
#     min = [eval_with(calib, domain[(k)][1]) for k in states]
#     max = [eval_with(calib, domain[(k)][2]) for k in states]
#     return Domain(states, min, max)
# end

# function get_grid(model::ASModel; options=Dict())
#     domain = get_domain(model)
#     d = length(domain.states)
#     options = get_options(model; options=Dict(:grid=>options))
#     grid_dict = options[:grid]
#     if :type in keys(grid_dict)
#         grid_dict[:tag] = grid_dict[:type]
#     end
#     if grid_dict[:tag] == :Cartesian
#         orders = get(grid_dict, :orders, [20 for i=1:d])
#         grid = CartesianGrid{d}(domain.min, domain.max, orders)
#         if length(orders)!=length(model.calibration[:states])
#             msg = string("Check the dimension of the matrix given in the yaml file, section: options-grid-orders. ",
#                          "Expected to be of dimension $([length(model.calibration[:states])])")
#             error(msg)
#         end
#     elseif grid_dict[:tag] == :Smolyak
#         mu = get(grid_dict, :mu, 3)
#         grid = SmolyakGrid{d}(domain.min, domain.max, mu)
#     elseif grid_dict[:tag] == :Random
#         n = get(grid_dict, :N, 200)
#         grid = RandomGrid{d}(domain.min, domain.max, n)
#     else
#         error("Unknown grid type.")
#     end
#     return grid
# end

# function get_exogenous_old(model::ASModel)
#     exo_dict = get(model.data,:exogenous,Dict{Symbol,Any}())
#     if length(exo_dict)==0
#         exo_dict = get(model.data[:options], :exogenous, Dict{Symbol,Any}())
#     end
#     exo = _symbol_dict(exo_dict)
#     calib = get_calibration(model)
#     exogenous = _build_exogenous_entry(exo, calib)
#     return exogenous
# end

# #


# function get_exogenous(model::AModel)
#     rdata = model.rdata
#     exo_dict = model.rdata[:exogenous]
#     cond = !(exo_dict.tag=="tag:yaml.org,2002:map")
#     if cond
#         # old style exogenous block
#         return get_exogenous_old(model)
#     else
#         syms = cat( [[Symbol(strip(e)) for e in split(k, ",")] for k in keys(exo_dict)]..., dims=1)
#         expected = model.symbols[:exogenous]
#         if (syms != expected)
#             msg = string("In 'exogenous' section, shocks must be declared in the same order as shock symbol. Found: ", syms, ". Expected: ", expected, ".")
#             throw(ErrorException(msg))
#         end
#         calibration = model.calibration.flat
#         processes = []
#         for k in keys(exo_dict)
#             v = exo_dict[k]
#             p = Dolang.eval_node(v, calibration, minilang, FromGreek())
#             push!(processes, p)
#         end
#         return ProductProcess(processes...)
#     end
# end


# function set_calibration!(model::ASModel, key::Symbol, value::Union{Real,Expr, Symbol})
#     model.data[:calibration][key] = value
# end


# mutable struct Model{ID} <: AModel{ID}

#     data::Dict{Symbol,Any}
#     rdata::YAML.Node
#     name::String           # weakly immutable
#     filename::String
#     symbols::OrderedDict{Symbol,Array{Symbol,1}}        # weakly immutable
#     equations      # weakly immutable
#     definitions    # weakly immutable
#     # functions      # weakly immutable
#     factories::Dict
#     calibration::ModelCalibration
#     exogenous
#     domain
#     grid
#     options

#     function Model{ID}(data::Dict{Symbol,Any}, rdata::YAML.Node; print_code::Bool=false, filename="<string>") where ID

#         model = new{ID}(data, rdata)
#         model.name = get_name(model)
#         model.filename = filename
#         model.symbols = get_symbols(model)
#         model.equations = get_equations(model)
#         model.definitions = get_definitions(model)
#         model.factories = Dict()
#         model.calibration = get_calibration(model)
#         model.exogenous = get_exogenous(model)
#         model.domain = get_domain(model)
#         model.options = get_options(model)
#         model.grid = get_grid(model)

#         factories = Dict{Symbol,Dolang.FlatFunctionFactory}()
#         for (eqtype,equations) in model.equations
#             spec = RECIPES[:dtcc][:specs][eqtype]
#             spec_eqs = spec[:eqs]
#             symbols = model.symbols
#             arguments = OrderedDict{Symbol,Union{Vector{Symbol},Vector{Tuple{Symbol,Int}}}}([
#                 Symbol(s)=>[(e,t) for e in symbols[Symbol(sg)]]
#                 for (sg,t,s) in spec_eqs
#             ])
#             arguments[:p] = [e[1] for e in arguments[:p]]
#             tdefs = get_definitions(model)

#             # if :target in keys(spec)
#             #     sg,t,s = spec[:target]
#             #     targets = [(e,t) for e in symbols[Symbol(sg)]]
#             #     peq = [Dolang.match_equality(eq) for eq in equations]
#             #     lhs = [eq[1] for eq in peq]
#             #     rhs = [eq[2] for eq in peq]
#             #     expressions = OrderedDict(zip(targets,rhs))
#             #     @assert Dolang.stringify.(lhs) == Dolang.stringify.(targets)
#             # else
#             #     expressions = OrderedDict( [Symbol("out_",i)=>eq_to_expr(eq) for (i,eq) in enumerate(equations)])
#             # end

#             # ff = Dolang.FlatFunctionFactory(expressions, arguments, tdefs; funname=eqtype )
#             # factories[eqtype] = ff
#             # code = Dolang.gen_generated_gufun(ff; dispatch=typeof(model))
#             # print_code && println("equation '", eqtype, "'", code)
#             # Core.eval(Dolo, code)
#         end

#         # model.factories = factories


#         # # TEMP: until we have a better method in Dolang
#         # # Create definitions
#         # vars = cat(model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls]; dims=1)
#         # args = OrderedDict(
#         #     :past => [(v,-1) for v in vars],
#         #     :present => [(v,0) for v in vars],
#         #     :future => [(v,1) for v in vars],
#         #     :params => model.symbols[:parameters]
#         # )
#         # fff = Dolang.FlatFunctionFactory(model.definitions, args, typeof(model.definitions)())
#         # code = Dolang.gen_generated_gufun(fff;dispatch=typeof(model), funname=:evaluate_definitions)
#         # Core.eval(Dolo,code)

#         return model
#     end

# end

# _numeric_mod_type(::Model{ID}) where {ID} = Model{ID}

# function Base.show(io::IO, model::Model)
#     print(io, "Model")
# end

# #### import functions

# function Model(url::AbstractString; print_code=false)
#     # it looks like it would be cool to use the super constructor ;-)
#     if match(r"(http|https):.*", url) != nothing
#         res = HTTP.request("GET", url)
#         txt = String(res.body)
#     else
#         txt = read(open(url), String)
#     end
#     data = _symbol_dict(load(txt, yaml_types))
#     rdata = Dolang.yaml_node_from_string(txt)
#     id = gensym()
#     fname = basename(url)
#     return Model{id}(data, rdata; print_code=print_code, filename=fname)
# end

# function set_calibration!(model::Model, key::Symbol, value::Union{Real,Expr, Symbol})
#     model.data[:calibration][key] = value
#     model.calibration = get_calibration(model)
#     model.exogenous = get_exogenous(model)
#     model.domain = get_domain(model)
#     # model.options = get_options(model) # we don't substitute calib here
#     model.grid = get_grid(model)
#     model.calibration
# end

# function set_calibration!(model::Model, values::AbstractDict{Symbol,T}) where T
#     for (key,value) in values
#         model.data[:calibration][key] = value
#     end
#     model.calibration = get_calibration(model)
#     model.exogenous = get_exogenous(model)
#     model.domain = get_domain(model)
#     # model.options = get_options(model) # we don't substitute calib here
#     model.grid = get_grid(model)
#     model.calibration
# end

# function set_calibration!(model::Model; kwargs...)
#     set_calibration!(model, Dict(kwargs))
# end

# """
# Imports the model from a yaml file specified by the `url` input
# parameter, and returns the corresponding `Dolo` model object.
# """
# function yaml_import(url; print_code::Bool=false)
#     Model(url; print_code=print_code)
# end


# function features(model::ASModel)
#     features = Dict{Symbol, Bool}()
#     features[:one_state] = (length(model.symbols[:states])==1)
#     features[:one_control] = (length(model.symbols[:controls])==1)
#     features[:one_dimensional] = features[:one_state] && features[:one_state]
#     it = Dolang.IncidenceTable( model.equations[:transition] )
#     features[:nonstochastic_transitions] = (length( intersect(it.by_date[0], model.symbols[:exogenous] ) )==0)
#     it = Dolang.IncidenceTable( [model.equations[:controls_lb], model.equations[:controls_ub]] )
#     features[:bounds_are_constant] = (length(it.by_var) == 0)
#     return features
# end
