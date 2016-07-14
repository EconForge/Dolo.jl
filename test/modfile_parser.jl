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


  # Remove lines beginning with comments
  filter!(x -> ~ismatch(r"^\s*%", x), lines)
  filter!(x -> ~ismatch(r"^\s*//", x), lines)

  # Remove lines with comments after text
  for ln = 1:length(lines)
    lines[ln] = split(lines[ln], r"%.*$")[1]
  end

  # Remove lines that are empty
  filter!(x -> (x != "\n"), lines)

  # Remove new lines
  for ln = 1:length(lines)
    lines[ln] = strip(lines[ln])
  end

  # Smoosh text back together
  text = join(lines)


  # Get variable names
  tmp = match(r"var (.*?);", text)
  variables = split(tmp[1], " ")

  # Get exogenous variable names
  tmp = match(r"varexo (.*?);", text)
  exogenous = split(tmp[1], " ")

  # Get parameter names
  tmp = match(r"parameters (.*?);", text)
  parameters = split(tmp[1], " ")

  # get equations, fill a list
  tmp  = match(r"model;(.*?)end;", text)
  equations = split(tmp[1], ";")
  equations = filter!(x -> length(x)>0, equations)


  # Get calibration values, fill a dictionary
  calibration = OrderedDict()
  tmp  = match(r"parameters(.*?);(.*)model", text)
  tmp = split(tmp[2], ";")
  tmp = filter!(x -> length(x)>0, tmp)
  for ln = 1:length(tmp)
    key = strip(match(r"(.*)=", tmp[ln])[1])
    entry = strip(match(r"=(.*)", tmp[ln])[1])
    calibration[key] =  entry
  end

  # Get initial values, fill a dictionary
  initval = OrderedDict()
  if contains(text, "initval")
    tmp  = match(r"initval;(.*?)end;", text)
    tmp = split(tmp[1], ";")
    tmp = filter!(x -> length(x)>0, tmp)
    for ln = 1:length(tmp)
      key = strip(match(r"(.*)=", tmp[ln])[1])
      entry = strip(match(r"=(.*)", tmp[ln])[1])
      initval[key] =  entry
    end
  end

  # Get end values, fill a dictionary
  endval = OrderedDict()
  if contains(text, "endval")
    tmp  = match(r"endval;(.*?)end;", text)
    tmp = split(tmp[1], ";")
    tmp = filter!(x -> length(x)>0, tmp)
    for ln = 1:length(tmp)
      key = strip(match(r"(.*)=", tmp[ln])[1])
      entry = strip(match(r"=(.*)", tmp[ln])[1])
      endval[key] =  entry
    end
  end

  # Get the calibrated shock values
  if contains(text, "shocks")
    tmp = match(r"shocks;(.*?)(.*?)end;", text)
    if contains(tmp[2], "stderr")
      shocks = fill("0", length(exogenous), length(exogenous))
      tmp = filter!(x-> (x != ""), split(tmp[2], ";"))
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


  close(f)

  return variables, exogenous, parameters, calibration, equations, initval, endval, shocks

end
