include("linter.jl")

filename = Pkg.dir("Dolo","linter/models/rbc_dtcc_iid.yaml")
filename


import DoloLinter
DL = DoloLinter


errors, warnings = DoloLinter.check_symbols(filename)



# import yaml structure:
import DoloLinter: Location, LinterException, get_loc, LinterWarning
node = DoloLinter.yaml_node_from_file(filename)

function check_name(d)
  if !("name" in keys(d)) || d["name"] == ""
    return false
  end
end

function check_model_type(d)
  if !("model_type" in keys(d)) || d["model_type"] == ""
    return false
  end
end

d= node

check_name(node)

check_model_type(d)
node[":type"]

loc = Location(node.start_mark, node.end_mark)
exc = LinterException("hi", loc, "<file>")


#### Some trials to find the repeating value in symbols ###~
# find
# - First occurence (location, symbol type)
# - Location of second or later occurence
# - Push as error if symbol != first_occurence_symbol

d = node

for (i, key)  in keys(node["symbols"])
  for (j, sd) in enumerate(node["symbols"]["key"].value)
    if count(c -> c == sd.value, collect(syms)) > 1
      floc = findfirst(syms, sd.value)
      if i > floc
        loc = get_loc(cat(1,values(d["symbols"])...)[floc])
        loc_reoccur = get_loc(sd)
        err_value = sd.value
        errtype = "Invalid symbol"
          msg = string(sd.value, "already declared")
          push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))
        end
    end
  end

end

d = node
syms = []
for sym in cat(1,values(d["symbols"])...)
  push!(syms, sym.value)
end

cat(1,values(d["symbols"])...)
d["symbols"]
cat(1,values(d["symbols"])...)

for (i, sd) in enumerate(cat(1,values(d["symbols"])...))
  if count(c -> c == sd.value, collect(syms)) > 1
    floc = findfirst(syms, sd.value)
    if i > floc
      loc_first = get_loc(cat(1,values(d["symbols"])...)[floc])
      loc = get_loc(sd)
      errvalue = sd.value
      errtype = "Invalid symbol"
      msg = string(sd.value, " already declared at line ", string(loc_first.start_mark.line))
      push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))
    end
  end
end

keys(d["symbols"])
cat(1,(d["symbols"]["parameters"].value)...)
values(d["symbols"])
values(d["symbols"])
cat(1,values(d["symbols"])...)
symnodef

for (i ,key) in enumerate(keys(d["symbols"]))
  for (j ,symnode) in enumerate(d["symbols"][key].value)
    if count(c -> c == symnode.value, collect(syms)) > 1
      floc = findfirst(syms, symnode.value)
      if symnode != cat(1,values(d["symbols"])...)[floc]
        loc = get_loc(symnode)
        errvalue = symnode.value
        errtype = "Invalid symbol"
        msg = string(symnode, " already declared in section ", )
        push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))
      end
    end
  end
end

function format_human(exc::Union{LinterException, LinterWarning})
    # Output should be a nice one line human-readable strings
    # colors, spacing, etc...
end

function format_human(errors::Vector{LinterException}, warnings::Vector{LinterWarning})
    # ouptputs all errors then all warnings, line-by-line
    # example from python's `dolo-lint --format=human examples/models/rbc_dtcc_iid.yaml`
    # error: 12, 30: Symbol 'beta' already declared as 'parameters'. (pos (12, 16))
end

function format_dict(exc::Union{LinterException, LinterWarning})
    # Output should be a dictionary like:
    # {"type": "error", "source": "examples/models/rbc_dtcc_iid.yaml", "range": [[11, 29], [11, 33]], "text": "Symbol 'beta' already declared as 'parameters'."}
end

function format_json(errors::Vector{LinterException}, warnings::Vector{LinterWarning})
    # Structure compatible with atom linter
    # see python's `dolo-lint --format=json examples/models/rbc_dtcc_iid.yaml` for an example
    # [{"type": "error", "source": "examples/models/rbc_dtcc_iid.yaml", "range": [[11, 29], [11, 33]], "text": "Symbol 'beta' already declared as 'parameters'. (pos (12, 16))"}]
end
