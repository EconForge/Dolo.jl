__precompile__(true)
module DoloLinter
    import YAML
    include("../src/linter.jl")
end

import DoloLinter: lint
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

lint(filename)
