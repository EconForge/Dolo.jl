using Dolang: sanitize
using Lerche: Tree, Token

using Dolo: ModelCalibration
import Dolang: solve_triangular_system

import YAML
import Dolang
import Dolang: yaml_node_from_string, yaml_node_from_file
using Dolang: SymExpr
import Dolang: LTree
import Dolang


using DataStructures: OrderedDict



# defind language understood in yaml files
language = minilang

mutable struct Model{ID, ExoT} <: AbstractModel{ExoT}

    data::YAML.MappingNode
    symbols::Dict{Symbol, Vector{Symbol}}
    calibration::ModelCalibration
    exogenous::ExoT
    domain::AbstractDomain
    factories
    definitions

    function Model{ID, ExoT}(data::YAML.Node) where ID where ExoT
        model = new{ID, ExoT}(data)
        model.calibration = get_calibration(model)
        model.symbols = get_symbols(model)
        model.exogenous = get_exogenous(model)
        model.domain = get_domain(model)
        
        factories = get_factories(model)
        model.factories = factories
        for (k,fact) in factories
            code = Dolang.gen_generated_gufun(fact; dispatch=typeof(model))
            # print_code && println("equation '", eq_type, "'", code)
            Core.eval(Dolo, code)
        end

        # Create definitions
        defs = get_definitions(model; stringify=true)
        model.definitions = defs
        definitions = OrderedDict{Symbol, SymExpr}( [(stringify(k),v) for (k,v) in defs])
        vars = cat(model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls]; dims=1)
        args = OrderedDict(
            :past => [Dolang.stringify(v,-1) for v in vars],
            :present => [Dolang.stringify(v,0) for v in vars],
            :future => [Dolang.stringify(v,1) for v in vars],
            :params => [Dolang.stringify(e) for e in model.symbols[:parameters]]
        )
        ff = Dolang.FunctionFactory(definitions, args, OrderedDict{Symbol, SymExpr}(), :definitions)
        code = Dolang.gen_generated_gufun(ff;dispatch=typeof(model), funname=:evaluate_definitions)
        Core.eval(Dolo,code)
        model.factories[:definitions] = ff
        

        return model
    end

end

id(::Model{ID}) where {ID} = ID


function check_exogenous_type(data)
    if !("exogenous" in keys(data))
        return EmptyProcess
    end
    exo = data["exogenous"]
    exo_types = [exo[k].tag for k in keys(exo)]

    iid_types = ("!UNormal", "!ULogNormal", "!Normal", "!MvNormal", "!ConstantProcess")
    mc_types = ("!MarkovChain", "!ConstantProcess")
    acorr_types = ("!AR1", "!VAR1", "!ConstantProcess")

    if exo_types ⊆ iid_types
        return   IIDExogenous
    elseif exo_types ⊆ mc_types
        return DiscreteProcess
    else
        return ContinuousProcess
    end
end


function Model(url::AbstractString; print_code=false)

    # it looks like it would be cool to use the super constructor ;-)
    if match(r"(http|https):.*", url) != nothing
        res = HTTP.request("GET", url)
        txt = String(res.body)
    else
        txt = read(open(url), String)
    end
    data = Dolang.yaml_node_from_string(txt)
    id = gensym()
    fname = basename(url)
    typ = check_exogenous_type(data)
    return Model{id, typ}(data)

end

yaml_import(filename::AbstractString) = Model(filename)

function Base.show(io::IO, model::Model)
    print(io, "Model")
end



function get_symbols(model::Model)
    syms = model.data[:symbols]
    symbols = OrderedDict{Symbol,Vector{Symbol}}()
    for k in keys(syms)
        symbols[Symbol(k)] = [Symbol(e.value) for e in syms[k]]
    end
    return symbols
end

function get_variables(model::Model)
    symbols = get_symbols(model)
    vars = cat(values(symbols)...; dims=1)
    dynvars = setdiff(vars, symbols[:parameters] )
    dynvars = union( dynvars, get_defined_variables(model) )
end

function get_defined_variables(model::Model)
    if !("definitions" in keys(model.data))
        return Symbol[]
    else
        # TODO: this is awfully redundant. This should be done "after" definitions are parsed
        try
            # OBSOLETE
            return [Symbol(e) for e in keys(model.data["definitions"])]
        catch
            defeqs = Dolang.parse_assignment_block(model.data["definitions"].value).children
            return [Symbol(eq.children[1].children[1].children[1].value) for eq in defeqs]

        end
    end
end

function get_name(model::Model)
    if "name" in keys(model.data)
        name = model.data[:name].value
    else
        name = "anonymous"
    end
    return name
end



function get_exogenous(model::AModel)

    rdata = model.data
    if !("exogenous" in keys(model.data))
        return nothing
    end
    exo_dict = model.data[:exogenous]
    cond = !(exo_dict.tag=="tag:yaml.org,2002:map")

    syms = cat( [[Symbol(strip(e)) for e in split(k, ",")] for k in keys(exo_dict)]..., dims=1)
    expected = model.symbols[:exogenous]
    if (syms != expected)
        msg = string("In 'exogenous' section, shocks must be declared in the same order as shock symbol. Found: ", syms, ". Expected: ", expected, ".")
        throw(ErrorException(msg))
    end
    calibration = model.calibration.flat
    processes = []
    for k in keys(exo_dict)
        v = exo_dict[k]
        p = Dolang.eval_node(v, calibration, minilang, ToGreek())
        push!(processes, p)
    end
    return ProductProcess(processes...)

end

function get_calibration(model::Model; kwargs...)
    calib = model.data[:calibration]
    # prep calib: parse to Expr, Symbol, or Number
    _calib  = OrderedDict{Symbol,Union{Expr,Symbol,Number}}()
    # add calibration for definitions
    for (k,v) in kwargs
        calib[k].value = string(v)
    end
    for k in keys(calib)
        x = Dolang.parse_equation(calib[k].value)
        e = Dolang.convert(Expr, x, stringify=false)
        _calib[Symbol(k)] = e
    end
    # so far _calib is a symbolic calibration
    calibration = solve_triangular_system(_calib)
    symbols = get_symbols(model)
    return ModelCalibration(calibration, symbols)
end

function set_calibration!(model::Model, key::Symbol, value)
    # TODO: set proper type for ScalarNode
    # this will fail is parameter wasn't defined before
    model.data[:calibration][key].value = string(value)
    model.calibration = get_calibration(model)
    model.exogenous = get_exogenous(model)
    model.domain = get_domain(model);
end

function set_calibration!(model::Model; values...)
    calib = model.data[:calibration]
    for (k,v) in values
        calib[k].value = string(v)
    end
    model.calibration = get_calibration(model)
    model.exogenous = get_exogenous(model)
    model.domain = get_domain(model);
end



function get_equation_block(model, eqname; stringify=true)
    variables = get_variables(model)
    variables_str = String[string(e) for e in variables]
    
    yml_node = model.data["equations"][eqname]

    if "tag:yaml.org,2002:str" == yml_node.tag
        block = Dolang.parse_equation_block(yml_node.value; variables=variables_str)
        equations = block.children
    else
        equations = [Dolang.parse_equation(line.value; variables=variables_str) for line in yml_node]
    end

    eqs = Union{Expr, Number, Symbol}[]
    eqs_lb = Union{Expr, Number, Symbol}[]
    eqs_ub = Union{Expr, Number, Symbol}[]    

    for eq in equations
        if eq.data == "double_complementarity"
            inequality = eq.children[2]
            # @assert inegality.data == "double_inequality"
            eq_lb = Dolang.convert(Expr, inequality.children[1]; stringify=stringify)
            eq_cc = Dolang.convert(Expr, inequality.children[2]; stringify=stringify)
            # TODO: check complementarities order
            eq_ub = Dolang.convert(Expr, inequality.children[3]; stringify=stringify)
            eq = eq.children[1]
            push!(eqs_lb, eq_lb)
            push!(eqs_ub, eq_ub)
        else
            push!(eqs_lb, :(-Inf) )
            push!(eqs_ub, :(Inf) )
        end
        if eq.data == "equality"
            eq = Tree("sub", [eq.children[2], eq.children[1]])          
        end
        # eq = Dolang.stringify(eq)
        push!(eqs, Dolang.convert(Expr, eq; stringify=stringify))
    end
    return (eqs, eqs_lb, eqs_ub)
end

function get_assignment_block(model, eqname; stringify=true)
    variables = get_variables(model)
    variables_str = String[string(e) for e in variables]
    
    yml_node = model.data["equations"][eqname]

    if "tag:yaml.org,2002:str" == yml_node.tag
        block = Dolang.parse_assignment_block(yml_node.value; variables=variables_str)
        equations = block.children
    else
        equations = [Dolang.parse_assignment(line.value; variables=variables_str) for line in yml_node]
    end


    # eqs = OrderedDict{Tuple{Symbol, Int64}, SymExpr}()
    eqs = OrderedDict()
    for eq in equations

        lhs = eq.children[1]
        rhs = eq.children[2]
        dest_s = Symbol(lhs.children[1].children[1].value)
        dest_d = parse(Int, lhs.children[2].children[1].value)
        
        dest = Dolang.stringify(dest_s, dest_d)
        
        # TODO : remove type instability
        if stringify
            dest = Dolang.stringify(dest_s, dest_d)
        else
            dest = (dest_s, dest_d)
        end
        
        expr = Dolang.convert(Expr, rhs; stringify=stringify)

        eqs[dest] = expr

    end
    return eqs
end


function get_domain(model::Model)::AbstractDomain

    states = model.symbols[:states]

    if !("domain" in keys(model.data))
        return EmptyDomain(states)
    end

    domain = model.data["domain"]
    
    kk = keys(domain)
    if "min" in kk
        # TODO: deprecation warning
        calib = model.calibration.flat
        min = [Dolang.eval_node(domain[(k)][1], calib) for k in states]
        max = [Dolang.eval_node(domain[(k)][2], calib) for k in states]
        return Domain(states, min, max)

    else
        calib = model.calibration.flat
        kk = tuple([Symbol(k) for k in keys(domain)]...)
        ss = tuple(states...)
        if kk != ss
            error("The declaration order in the domain section ($kk) must match the symbols orders ($ss)")
        end
        vals = [Dolang.eval_node(domain[k], calib) for k in kk]
        min = [e[1] for e in vals]
        max = [e[2] for e in vals]
        return Domain(states, min, max)
    end

end

function get_definitions(model::Model; tshift=0, stringify=false) # ::OrderedDict{Tuple{Symbol,Int}}
    
    # return OrderedDict{Tuple{Symbol,Int64},Union{Expr, Number, Symbol}}()

    # parse defs so values are Expr
    if !("definitions" in keys(model.data))
        return OrderedDict{Tuple{Symbol,Int64},Union{Expr, Number, Symbol}}()
    end

    if model.data["definitions"].tag == "tag:yaml.org,2002:map"

        # LEGACY
        defs = Dict()

        dynvars = string.( cat(get_variables(model), keys(defs)...; dims=1) )

        _defs =  OrderedDict{Tuple{Symbol,Int64},Union{Expr, Number, Symbol}}()  # {Symbol,Int}()

        for (kkk,v) in defs
            kk = Dolang.parse_equation(kkk)
            if kk.data == "symbol"
                k = kk.children[1].value
            elseif kk.data=="variable"
                k = kk.children[1].children[1].value
            else
                error("Invalid key.")
            end
            vv = Dolang.parse_equation(v.value; variables=(dynvars))::LTree
            if tshift != 0
                vv = Dolang.time_shift(vv, tshift)::LTree
            end
            _defs[(Symbol(k),tshift)] = Dolang.convert(Expr, vv; stringify=stringify)::SymExpr
        end

        return _defs

    else

        @assert model.data["definitions"].tag == "tag:yaml.org,2002:str"
        _defs =  OrderedDict{Tuple{Symbol,Int64},Union{Expr, Number, Symbol}}()  # {Symbol,Int}()
        equations = Dolang.parse_assignment_block(model.data["definitions"].value).children
        # TODO : check that equations are correct
        for eq in equations
            dest, val = eq.children
            # TODO : create conveniences functions in dolang to avoid the follwowing
            name = dest.children[1].children[1].value
            # TODO : check that name is already defined at this stage
            t = parse(Int, dest.children[2].children[1].value)

            @assert t == 0

            vv = Dolang.time_shift(val, tshift)

            _defs[(Symbol(name),tshift)] = Dolang.convert(Expr, vv; stringify=stringify)::SymExpr
        end
        return _defs
    end

end


function get_options(model::Model)
    cc = model.calibration.flat
    if "options" in keys(model.data)
        return Dolang.eval_node(model.data["options"], cc, language)
    else
        return Dict()
    end
end

import Dolang: FunctionFactory
import Dolang: stringify

function get_factories(model::Model)

    facts = Dict()
    for eq_type in ["transition", "arbitrage"]
        if eq_type == "transition"
            factories = (get_factory(model, eq_type), )
        else
            factories = get_factory(model, eq_type)
        end
        for f in factories
            facts[f.funname] = f
        end
    end

    return facts

end

function get_factory(model::Model, eq_type::String)
    if eq_type == "arbitrage"
        defs_0 = get_definitions(model; stringify=true)
        defs_1 = get_definitions(model; tshift=1, stringify=true)
        definitions = OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in merge(defs_0, defs_1)])

        eqs, eq_lb, eq_ub = get_equation_block(model, eq_type)

        symbols = get_symbols(model)

        #TODO : fix crazy bug: it doesn't work without the trailing underscore !
        equations = OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
            Symbol(string("out_", i, "_")) => eqs[i]
            for i =1:length(eqs)
        )
        arguments = OrderedDict(
            :m => [stringify(e,0) for e in symbols[:exogenous]],
            :s => [stringify(e,0) for e in symbols[:states]],
            :x => [stringify(e,0) for e in symbols[:controls]],
            :M => [stringify(e,1) for e in symbols[:exogenous]],
            :S => [stringify(e,1) for e in symbols[:states]],
            :X => [stringify(e,1) for e in symbols[:controls]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )
        
        ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))

        equations_lb = OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
                Symbol(string("out_", i, "_")) => eq_lb[i]
                for i =1:length(eqs)
            )
        equations_ub = OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
                Symbol(string("out_", i,"_")) => eq_ub[i]
                for i =1:length(eqs)
            )
        arguments_bounds = OrderedDict(
            :m => [stringify(e,0) for e in symbols[:exogenous]],
            :s => [stringify(e,0) for e in symbols[:states]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )

        # definitions: we should remove definitions depending on controls
        # definitions = OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in defs_0] )
        definitions = Dict{Symbol, SymExpr}()
        ff_lb = FunctionFactory(equations_lb, arguments_bounds, definitions, :controls_lb)
        ff_ub = FunctionFactory(equations_ub, arguments_bounds, definitions, :controls_ub)

        return (ff, ff_lb, ff_ub)
           
    else
        # defs_0  = get_definitions(model; stringify=true)
        defs_m1 = get_definitions(model; tshift=-1, stringify=true)

        definitions = OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in  defs_m1 ])

        equations = get_assignment_block(model, eq_type)
        symbols = get_symbols(model)
        arguments = OrderedDict(
            :m => [stringify(e,-1) for e in symbols[:exogenous]],
            :s => [stringify(e,-1) for e in symbols[:states]],
            :x => [stringify(e,-1) for e in symbols[:controls]],
            :M => [stringify(e,0) for e in symbols[:exogenous]],
            :p => [stringify(e) for e in symbols[:parameters]]
        )
        ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
        return ff
    end

end


function get_options(model::AModel; options=Dict())
    if "options" in keys(model.data)
        return model.data["options"]
        # opts = model.data["options"]
        # return rmerge(opts, options)
    else
        return Dict()
    end
end

using Dolo: CartesianGrid

function get_discretization_options(model::AModel)

    domain = get_domain(model)
    d = length(domain.states)

    c = model.calibration.flat

    if ("options" in keys(model.data)) && ( "discretization" in keys(model.data["options"]) )
        gopt = model.data["options"]["discretization"]
        return Dolang.eval_node(gopt, c)
        

    else
        return Dict(
            :endo=>Dict(),
            :exo=>Dict()
        )
    end


        
end


function discretize(model; kwargs...)

    opts = get_discretization_options(model; kwargs...)
    
    opts_endo = merge(opts[:endo], get(kwargs, :endo, Dict()) )
    opts_exo = merge(opts[:exo], get(kwargs, :exo, Dict()) )

    grid_endo = Dolo.discretize(model.domain;  opts_endo...) 
    dprocess = Dolo.discretize(model.exogenous;  opts_exo...) 
    grid_exo = dprocess.grid
    grid = ProductGrid(grid_exo, grid_endo)
    return grid, dprocess


    # dprocess = discretize(model.exogenous)
#     return ProductGrid(dprocess.grid, grid), dprocess
    # exo_grid = dprocess.
    # return ProductGrid(grid, exo_grid)
    # return dprocess.grid, grid, dprocess
end

