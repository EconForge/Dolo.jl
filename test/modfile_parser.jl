#=
This file reads Dynare .mod files, and extracts the following information:
  variables
  exo_variables
  parameters
  calibration
  model
  initval
  endval
  shocks
=#

using DataStructures

function modfile_parser(mod_file_name)
  f = open(mod_file_name)
  lines = readlines(f)  # reads in line by line
  close(f)

  # Remove lines beginning with comments (%, //, /*, \*)
  filter!(x -> ~ismatch(r"^\s*%", x), lines)
  filter!(x -> ~ismatch(r"^\s*\/\/", x), lines)
  filter!(x -> ~ismatch(r"^\s*\/\*", x), lines)


  # Remove end-of-line comments
  for ln = 1:length(lines)
    lines[ln] = split(lines[ln], r"%.*$")[1]    # %  comments
    lines[ln] = split(lines[ln], r"\/\/.*$")[1]  # \\ comments
    lines[ln] = split(lines[ln], r"/\*.*$")[1]   # \* comments
  end


  # Smoosh text back together
  text = join(lines)

  # Remove new lines, spaces, carriage returns,
  text = replace(text, "\t", " ")
  text = replace(text, "\r", " ")
  text = replace(text, "\n", " ")

  # Get variable names
  tmp = match(r"var\s(.*?);", text)
  variables = split(tmp[1], " ")
  filter!(x -> ismatch(r"\S", x), variables) # Remove lines if no chars present

  # Get exogenous variable names
  tmp = match(r"varexo\s(.*?);", text)
  shocks = split(tmp[1], " ")
  filter!(x -> ismatch(r"\S", x), shocks) # Remove lines if no chars present

  # Get parameter names
  tmp = match(r"parameters\s(.*?);", text)
  parameters = split(tmp[1], " ")
  filter!(x -> ismatch(r"\S", x), parameters) # Remove lines if no chars present

  # get equations, fill a list
  tmp  = match(r"model;\s(.*?)end;", text)
  equations = split(tmp[1], ";")
  filter!(x -> ismatch(r"\S", x), equations) # Remove lines if no chars present
  for ln = 1:length(equations)
    equations[ln] = strip(equations[ln])
  end

  # Get calibration values, fill a dictionary
  parameter_values = OrderedDict()
  tmp  = match(r"parameters\s(.*?);(.*)model", text)
  tmp = split(tmp[2], ";")
  filter!(x -> ismatch(r"\S", x), tmp) # Remove lines if no chars present
  for ln = 1:length(tmp)
    key = strip(match(r"(.*)=", tmp[ln])[1])
    entry = strip(match(r"=(.*)", tmp[ln])[1])
    parameter_values[key] =  entry
  end

  # Get initial values, fill a dictionary
  initval = OrderedDict()
  tmp  = match(r"initval;(.*?)end;", text)
  if tmp != nothing
    tmp = split(tmp[1], ";")
    filter!(x -> ismatch(r"\S", x), tmp) # Remove lines if no chars present
    for ln = 1:length(tmp)
      key = strip(match(r"(.*)=", tmp[ln])[1])
      entry = strip(match(r"=(.*)", tmp[ln])[1])
      initval[key] =  entry
    end
  end



  # Get end values, fill a dictionary
  endval = OrderedDict()
  tmp  = match(r"endval;(.*?)end;", text)
  if tmp != nothing
    tmp = split(tmp[1], ";")
    filter!(x -> ismatch(r"\S", x), tmp) # Remove lines if no chars present
    for ln = 1:length(tmp)
      key = strip(match(r"(.*)=", tmp[ln])[1])
      entry = strip(match(r"=(.*)", tmp[ln])[1])
      endval[key] =  entry
    end
  end



  # Get the calibrated shock values, fill matrix
  shockvaldict = Dict()
  shock_matrix = fill("0", length(shocks), length(shocks))
  tmp = match(r"shocks;(.*?)(.*?)end;", text)
  if tmp != nothing
    if contains(tmp[2], "stderr")
      tmpkey = matchall(r"var(.*?);", tmp[2])
      tmpentry = matchall(r"stderr(.*?);", tmp[2])
      for ln = 1:length(tmpkey)
       key = strip(match(r"var\s(.*)", tmpkey[ln])[1])
       key = match(r"(.*);", key)[1]
       entry = strip(match(r"stderr\s(.*)", tmpentry[ln])[1])
       entry = match(r"(.*);", entry)[1]
       shockvaldict[key] =  entry
      end
    else
      tmp = split(tmp[2], ";")
      filter!(x -> ismatch(r"\S", x), tmp) # Remove lines if no chars present
      for ln = 1:length(tmp)
        key = match(r"var\s(.*)\s=", tmp[ln])[1]
        entry = match(r"=\s(.*)", tmp[ln])[1]
        shockvaldict[key] =  entry
      end
    end
    # Fill matrix in order of shocks in "shocks" dictionary
    for ln = 1:length(shocks)
      if haskey(shockvaldict, shocks[ln])
        shock_matrix[ln, ln] = shockvaldict[shocks[ln]]
      else
        shock_matrix[ln, ln] = "0.0"
      end
    end
  end

  data = Dict()
  data["symbols"] = Dict("variables"=>variables, "shocks"=>shocks, "parameters"=>parameters)
  data["equations"] = equations
  data["calibration"] = merge(parameter_values, initval)
  distribution = Dict("tag"=>"Normal", "sigma"=>shock_matrix)
  data["options"] = Dict("distribution"=>distribution)
  data["optional"] = Dict("endval"=>endval)

  return data

end
