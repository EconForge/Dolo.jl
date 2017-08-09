include("linter.jl")

filename = Pkg.dir("Dolo","linter/models/rbc_dtcc_iid.yaml")
filename


import DoloLinter
DL = DoloLinter

model_symbol_types = keys(d[:symbols])
errors, warnings = DoloLinter.check_symbols(filename)



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
