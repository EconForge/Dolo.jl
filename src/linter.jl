function yaml_node_from_string(fn::AbstractString)

    # Didn't find how to access the top node more easily.
    #
    yml_types = Dict{String,Function}()
    yml_types["!Cartesian"] = ((c,n)->n)
    yml_types["!Smolyak"] = ((c,n)->n)
    yml_types["!Normal"] = ((c,n)->n)
    yml_types["!MarkovChain"] = ((c,n)->n)
    yml_types["!Product"] = ((c,n)->n)
    yml_types["!PoissonProcess"] = ((c,n)->n)
    yml_types["!DeathProcess"] = ((c,n)->n)
    yml_types["!AgingProcess"] = ((c,n)->n)
    yml_types["!VAR1"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:null"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:bool"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:int"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:float"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:binary"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:timestamp"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:omap"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:pairs"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:set"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:str"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:seq"] = ((c,n)->n)
    yml_types["tag:yaml.org,2002:map"] = ((c,n)->n)


    return YAML.load(fn, yml_types)
end


function yaml_node_from_file(fn::AbstractString)
    txt = open(f->read(f,String), fn)
    txt = replace(txt, "\r"=>"")
    return yaml_node_from_string(txt)
end

mutable struct Location
    start_mark::YAML.Mark
    end_mark::YAML.Mark
end

mutable struct MissingElement <: Exception
    path::Vector{AbstractString}
    k::Integer # index of first missing element
    loc::Location
end

function Base.getindex(d::YAML.MappingNode, l::AbstractString...)
    el = d
    for k=1:length(l)
        # if !(typeof(el)<:YAML.MappingNode)
        #     println("This is not going to fly")
        if !(l[k] in keys(el))
            loc = Location(el.start_mark, el.end_mark)
            throw(MissingElement(collect(l), k, loc))
        else
            el = el[l[k]]
        end
    end
    return el
end
function Base.getindex(d::YAML.MappingNode, l::Symbol...)
    return Base.getindex(d, [string(e) for e in l]...)
end


mutable struct LinterException <: Exception
    msg::AbstractString
    loc::Location
    src::AbstractString
end

mutable struct LinterWarning <: Exception
    errvalue::AbstractString
    errtype::AbstractString
    msg::AbstractString
    loc::Location
    src::AbstractString
end

#LinterException(msg::AbstractString, loc::Location) = LinterException(msg, loc, "<string>")

LinterWarning(msg::AbstractString, loc::Location) = LinterWarning(msg, loc, "<string>", errtype)

function get_loc(d)
    Location(d.start_mark, d.end_mark)
end

function check_name(d)
  if !("name" in keys(d)) || typeof(d["name"].value ) != String || d["name"].value == ""
    return false
  end
end


function check_symbol_validity(sym)
  # Not sure if it is necesarry to check if that's a string: yaml structure saves by default as string(??)
  if  typeof(sym) != String || match(r"^[a-zA-Z0-9_]*$",sym) == nothing || string(sym)[1:1] == "_"
    return false
  end
end

#function check_equations(d::YAML.MappingNode, filename="<string>")
#Known_equation_types -> ["transition", "arbitrage", "value", "felicity", "expectation"]
#end

function check_model_sections(d::YAML.MappingNode, filename="<string>")
  errors = LinterWarning[]
  warnings = LinterWarning[]

  # Can add other reauired models sections: symbols, equations ...
  if check_name(d) == false
    if !("name" in keys(d))
      msg = ""
      errvalue = "name"
      errtype = "Missing model part"
      loc = get_loc(d)
    elseif typeof(d["name"].value ) != String || d["name"].value == ""
      errvalue = string(d["name"].value)
      errtype = "Invalid model name"
      loc = get_loc(d["name"])
      msg = "Model name should be non-empty string"
          # push!(errors/warnings, LinterWarning(errvalue, errtype, msg, loc, src)
    end

    push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))

  end

  return [errors, warnings]
end

function check_calibration(d::YAML.MappingNode, filename="<string>")
  errors = LinterWarning[]
  warnings = LinterWarning[]

  model_symbol_types = keys(d[:symbols])
  model_symbols = []
  model_symbols_vec = []
  for key in keys(d[:symbols])
      n = length(d[:symbols][key].value)
      for i=1:n
          # would be cooler to do directly for sym in syms[vg]
          sym = d[:symbols][key].value[i]
          push!(model_symbols_vec, (key, sym.value, sym))
          push!(model_symbols, sym.value)
      end
  end
  model_symbols = cat(model_symbols, keys(d[:definitions]); dims=1)


  for (i, key) in enumerate(keys(d["calibration"]))
    if !(key in model_symbols)

      errvalue = key
      errtype = "Calibration not in symbols"
      loc = get_loc(d["calibration"].value[i][1])
      msg = string(key, " is not declared in the symbols section ")
      push!(warnings, LinterWarning(errvalue, errtype, msg, loc, filename))

    elseif count(c -> c == key, collect( keys(d["calibration"]))) > 1
      floc = something(findfirst(isequal(key), keys(d["calibration"])), 0)
      loc_first = get_loc(d["calibration"].value[floc][1])
      if i > floc
        errvalue = key
        errtype = "Invalid calibration"
        loc = get_loc(d["calibration"].value[i][1])
        msg = string(key, " already declared in calibration section at line ",loc_first.start_mark.line, ", column ",loc_first.start_mark.column)
        push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))

      end
    end



  end

  return [errors, warnings]
end

function check_equations(d::YAML.MappingNode, filename="<string>")
  errors = LinterWarning[]
  warnings = LinterWarning[]

  # check symbol names:
  required_equation_types = ["transition"]
  optional_equation_types = ["arbitrage", "value", "felicity", "expectation", "direct_response"]
  known_equation_types = cat(required_equation_types, optional_equation_types; dims=1)

  model_equation_types = keys(d[:equations])

  node_transition = YAML.ScalarNode[]
  node_arbitrage = YAML.ScalarNode[]
  for eq in d[:equations].value
    if eq[1].value == "transition"
      node_transition = eq[1]
    elseif eq[1].value == "arbitrage"
      node_arbitrage = eq[1]
    end

  end



  for s in required_equation_types
      if !(s in model_equation_types)
          errvalue = string(s)
          errtype = "Missing equation type"
          msg = ""
          loc = get_loc(d[:equations])
          # push!(errors/warnings, LinterWarning(errvalue, errtype, msg, loc, src)
          push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))
      end
  end

  for (i,sg) in enumerate(keys(d["equations"]))
      if !(sg in cat(known_equation_types; dims=1))
          errvalue = string(sg)
          errtype = "Unknown equation type"
          msg = ""
          # get precise location
          loc = get_loc(d["equations"].value[i][1])
          push!(warnings, LinterWarning(errvalue, errtype, msg, loc, filename))
      end
  end


  if ("transition" in model_equation_types) && ("states" in keys(d[:symbols]))
    n_transition = length(d[:equations]["transition"].value)
    n_states = length(d[:symbols]["states"].value)



    if n_transition != n_states

      errvalue = string(node_transition.value)
      errtype = "Invalid number of equations"
      msg = "Number of transition equations should be equal to number of states"
      # get precise location
      loc = get_loc(node_transition)

      push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))
    end




  end


  if ("arbitrage" in model_equation_types) && ("controls" in keys(d[:symbols]))
    n_arbitrage = length(d[:equations]["arbitrage"].value)
    n_controls = length(d[:symbols]["controls"].value)


    if n_arbitrage != n_controls
      errvalue = string(node_arbitrage.value)
      errtype = "Invalid number of equations"
      msg = "Number of arbitrage equations should be equal to number of controls"
      # get precise location
      loc = get_loc(node_arbitrage.value)
      push!(warnings, LinterWarning(errvalue, errtype, msg, loc, filename))
    end


  end

  return [errors, warnings]
end



function check_symbols(d::YAML.MappingNode, filename="<string>")

    ### so far, this is just an example

    errors = LinterWarning[]
    warnings = LinterWarning[]


    # check symbol names:
    required_symbol_types = ["states", "controls", "exogenous", "parameters"]
    optional_symbol_types = [ "values", "rewards", "expectations"]
    known_symbol_types = cat(required_symbol_types, optional_symbol_types; dims=1)

    model_symbol_types = keys(d[:symbols])
    model_symbols = []
    model_symbols_vec = []
    for key in keys(d[:symbols])
        n = length(d[:symbols][key].value)
        for i=1:n
            # would be cooler to do directly for sym in syms[vg]
            sym = d[:symbols][key].value[i]
            push!(model_symbols_vec, (key, sym.value, sym))
            push!(model_symbols, sym.value)
        end
    end




    for s in required_symbol_types
        if !(s in model_symbol_types)
            errvalue = string(s)
            errtype = "Missing symbol type"
            msg = ""
            loc = get_loc(d[:symbols])
            # push!(errors/warnings, LinterWarning(errvalue, errtype, msg, loc, src)
            push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))
        end
    end

    for (i,sg) in enumerate(keys(d["symbols"]))
        if !(sg in cat(known_symbol_types; dims=1))
            errvalue = string(sg)
            errtype = "Unknown symbol type"
            msg = ""
            # get precise location
            loc = get_loc(d["symbols"].value[i][1])
            push!(warnings, LinterWarning(errvalue, errtype, msg, loc, filename))
        end
    end

    for (i,m) in enumerate(values(d["symbols"]))
      for (j, n) in enumerate(values(d["symbols"])[i])
        sym =n.value
        if check_symbol_validity(sym) == false
          errvalue = string(sym)
          errtype = "Invalid symbol"
          loc = get_loc(n)
          # get precise location
          if  typeof(sym) != String
            msg = "symbol should be a string"
          elseif match(r"^[a-zA-Z0-9_]*$",sym) == nothing
            msg = "symbol should be an 'alphanumeric' string"
          elseif string(sym)[1:1] == "_"
            msg = "symbol should not start with a number or an underscore"
          end
          push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))
        end
      end
    end

    for (i,symnode) in enumerate(model_symbols_vec)
      if count(c -> c == symnode[2], collect(model_symbols)) > 1
        floc = something(findfirst(isequal(symnode[2]), model_symbols), 0)
        loc_first = get_loc(model_symbols_vec[floc][3])
        if i > floc
          errvalue = symnode[2]
          errtype = "Invalid symbol"
          loc = get_loc(symnode[3])
          msg = string(symnode[2], " already declared as '", symnode[1] , " at line ",loc_first.start_mark.line, ", column ",loc_first.start_mark.column)
          push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))
        end
      end
    end

    for (i, sym) in enumerate(keys(d["definitions"]))

      if check_symbol_validity(sym) == false
        errvalue = string(sym)
        errtype = "Invalid definition"
        loc = get_loc(d["definitions"].value[i][1])
        # get precise location
        if  typeof(sym) != String
          msg = "symbol should be a string"
        elseif match(r"^[a-zA-Z0-9_]*$",sym) == nothing
          msg = "symbol should be an 'alphanumeric' string"
        elseif string(sym)[1:1] == "_"
          msg = "symbol should not start with a number or an underscore"
        end
        push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))
      end

    end

    for (i, key) in enumerate(keys(d["definitions"]))
      if (key in model_symbols)
        floc = something(findfirst(isequal(key), model_symbols), 0)
        loc_first = get_loc(model_symbols_vec[floc][3])
        declared_sym_type = model_symbols_vec[floc][1]

        errvalue = key
        errtype = "Invalid definition"
        loc = get_loc(d["definitions"].value[i][1])
        msg = string(key, "  already declared in symbols section as '", declared_sym_type , "' at line",loc_first.start_mark.line, ", column ",loc_first.start_mark.column)
        push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))

      elseif count(c -> c == key, collect( keys(d["definitions"]))) > 1
        floc = something(findfirst(isequal(key), keys(d["definitions"])), 0)
        loc_first = get_loc(d["definitions"].value[floc][1])
        if i > floc
          errvalue = key
          errtype = "Invalid definition"
          loc = get_loc(d["definitions"].value[i][1])
          msg = string(key, " already declared in definitions section at line ",loc_first.start_mark.line, ", column ",loc_first.start_mark.column)
          push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))

        end
      end
    end

    for (i,symnode) in enumerate(model_symbols_vec)
      if !(symnode[2] in keys(d["calibration"]))
          errvalue = symnode[2]
          errtype = "Symbol not calibrated"
          loc = get_loc(symnode[3])
          msg = string(symnode[2], " not found in calibration section")
          push!(warnings, LinterWarning(errvalue, errtype, msg, loc, filename))
      end
    end

    return [errors, warnings]

end




function check_equations(fn::AbstractString)
    d = yaml_node_from_file(fn)
    errs, wars = check_equations(d, fn)
end

function check_calibration(fn::AbstractString)
    d = yaml_node_from_file(fn)
    errs, wars = check_calibration(d, fn)
end


function check_symbols(fn::AbstractString)
    d = yaml_node_from_file(fn)
    errs, wars = check_symbols(d, fn)
end

function check_model_sections(fn::AbstractString)
    d = yaml_node_from_file(fn)
    errs, wars = check_model_sections(d, fn)
end

function check_model(fn::AbstractString)
    d = yaml_node_from_file(fn)
    errs, wars = check_model(d, fn)
end

function print_error(err::LinterWarning)
printstyled("error "; color=:light_red)
print("at line ",err.loc.start_mark.line, ", column ",err.loc.start_mark.column, " : ")
print(err.errtype)
if err.msg != ""
  printstyled(" '",err.errvalue,"' "; color=:light_green)
  println(">> ", err.msg)
else
  printstyled(" '",err.errvalue,"' "; color=:light_green)
  println()
end

end

function print_warning(err::LinterWarning)
printstyled("warning "; color=:light_blue)
print("at line ",err.loc.start_mark.line, ", colummn ",err.loc.start_mark.column, " : ")
print(err.errtype)
if err.msg != ""
  printstyled(" '",err.errvalue,"' "; color=:light_green)
  println(">> ", err.msg)
else
  printstyled(" '",err.errvalue,"' "; color=:light_green)
  println()
end
end


function format_human(errors::Vector{LinterWarning}, warnings::Vector{LinterWarning})
for err in cat(errors; dims=1)
    print_error(err)
end

for err in cat(warnings; dims=1)
    print_warning(err)
end
end



function lint(filename::AbstractString; format=:human)

  funnames = [check_model_sections, check_symbols, check_equations, check_calibration]
  errors = LinterWarning[]
  warnings = LinterWarning[]

  for fun in funnames
    errors2, warnings2 = fun(filename)
    errors = cat(errors, errors2; dims=1)
    warnings = cat(warnings, warnings2; dims=1)
  end

  if format == :human
      format_human(errors,warnings)
  else
      println("Format '", format, "' not implemented (yet).")
  end
  return cat(errors, warnings; dims=1)
end
