include("linter.jl")
import DoloLinter: check, format_human


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

errors, warnings = check(filename)

if format == "human"
    return format_human(errors,warnings)
else
    println("Format '", format, "' not implemented (yet).")
end
