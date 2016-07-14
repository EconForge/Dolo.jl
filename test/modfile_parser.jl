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
  tmp = match(r"var (.*?);", text)
  variables = split(tmp[1], " ")
  filter!(x -> ismatch(r"\S", x), variables) # Remove lines if no chars present

  # Get exogenous variable names
  tmp = match(r"varexo (.*?);", text)
  exogenous = split(tmp[1], " ")
  filter!(x -> ismatch(r"\S", x), exogenous) # Remove lines if no chars present

  # Get parameter names
  tmp = match(r"parameters (.*?);", text)
  parameters = split(tmp[1], " ")
  filter!(x -> ismatch(r"\S", x), parameters) # Remove lines if no chars present

  # get equations, fill a list
  tmp  = match(r"model;(.*?)end;", text)
  equations = split(tmp[1], ";")
  filter!(x -> ismatch(r"\S", x), equations) # Remove lines if no chars present
  for ln = 1:length(equations)
    equations[ln] = strip(equations[ln])
  end

  # Get calibration values, fill a dictionary
  calibration = OrderedDict()
  tmp  = match(r"parameters(.*?);(.*)model", text)
  tmp = split(tmp[2], ";")
  filter!(x -> ismatch(r"\S", x), tmp) # Remove lines if no chars present
  for ln = 1:length(tmp)
    key = strip(match(r"(.*)=", tmp[ln])[1])
    entry = strip(match(r"=(.*)", tmp[ln])[1])
    calibration[key] =  entry
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

  # Get the calibrated shock values
  tmp = match(r"shocks;(.*?)(.*?)end;", text)
  if tmp != nothing
    if contains(tmp[2], "stderr")
      shocks = fill("0", length(exogenous), length(exogenous))
      tmp = split(tmp[2], ";")
      filter!(x -> ismatch(r"\S", x), tmp) # Remove lines if no chars present
      # tmp = filter!(x-> (x != ""), split(tmp[2], ";"))
      cnt = 1
      for ln = 1:length(tmp)
        if contains(tmp[ln], "stderr")
          shock_val = match(r"stderr (.*)", tmp[ln])[1]
          shocks[cnt, cnt] = shock_val
          cnt += 1
        end
      end
    else
      shocks = []
      tmp = filter!(x-> (x != ""), split(tmp[2], ";"))
      for ln = 1:length(tmp)
        shock_val = match(r"=(.*)", tmp[ln])[1]
        push!(shocks, shock_val)
      end
    end
  end

  return variables, exogenous, parameters, calibration, equations, initval, endval, shocks

end
