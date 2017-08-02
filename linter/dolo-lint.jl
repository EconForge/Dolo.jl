include("linter.jl")
import DoloLinter: check


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
    for err in cat(1,errors,warnings)
        println(err)
    end
else
    throw(Exception("Format '", format, "' not implemented (yet)."))
end
