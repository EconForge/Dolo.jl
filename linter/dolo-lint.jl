__precompile__(true)
module DoloLinter
    import YAML
    include("../src/linter.jl")
end

import DoloLinter: check_model, check_symbols, format_human, check_model_sections, check_equations, check_calibration

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--format"
        arg_type = String
        default = "human"
        help = "output format"
    "filename"
        help = "filename"
        required = true
end

parsed_args = parse_args(ARGS, s)

filename = parsed_args["filename"]
format = parsed_args["format"]

errors, warnings = check_model_sections(filename)

errors2, warnings2 = check_symbols(filename)
append!(errors, errors2)
append!(warnings, warnings2)

errors2, warnings2 = check_equations(filename)
append!(errors, errors2)
append!(warnings, warnings2)

errors2, warnings2 = check_calibration(filename)
append!(errors, errors2)
append!(warnings, warnings2)


if format == "human"
    return format_human(errors,warnings)
else
    println("Format '", format, "' not implemented (yet).")
end
