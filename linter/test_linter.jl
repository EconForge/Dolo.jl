import YAML
import Dolang
path = Pkg.dir("Dolo","src","linter.jl")
include(path)



filename = Pkg.dir("Dolo","linter/models/rbc_dtcc_iid.yaml")
filename

#DL = DoloLinter
d = yaml_node_from_file(filename)
errswarrs = lint(filename)


model_equation_types = keys(d[:equations])
model_equations = []
model_equations_vec = []
for key in keys(d[:equations])
    n = length(d[:equations][key].value)
    for i=1:n
        # would be cooler to do directly for sym in syms[vg]
        eq = d[:equations][key].value[i]
        push!(model_equations_vec, (key, eq.value, eq))
        push!(model_equations, eq.value)
    end
end


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


for eq in model_equations_vec
  if eq[1] == "arbitrage"

    equ , comp = split(eq[2], " | ")
    for sym in keys(Dolang.IncidenceTable(parse(equ)).by_var)

    end
  end
end

import Dolo: RECIPES


errors = LinterWarning[]
warnings = LinterWarning[]


eqg = "arbitrage"
spec = RECIPES[:dtcc][:specs][:arbitrage][:eqs]
Symbol_set = Dict(eqg => Set(([el[2] for el in spec if el[1]==eqg])) for eqg in Set(el[1] for el in spec))

for eq_ in d[:equations][eqg].value
      eqs = eq_.value
      equ , comp = split(eqs, " | ")
      inc = Dolang.IncidenceTable(parse(equ)).by_var
      for sym in keys(inc)
        sym = string(sym)
        if !(sym in model_symbols)  && !(string(sym) in keys(d[:definitions]))
          errvalue = sym
          errtype = "Unkown Symbol"
          loc = get_loc(eq_)
          msg = string(sym, "  in the '", eqg , "' equation is not defined in symbols or definitions sections")
          println(sym)
          push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))

        elseif  (string(sym) in model_symbols)
          floc = findfirst(model_symbols, sym)
          declared_sym_type = model_symbols_vec[floc][1]
          println(declared_sym_type)

          if !(inc[sym] in Symbol_set[declared_sym_type])
            errvalue = sym
          ##  errtype = "Timing convention error"
          ##  loc = get_loc(eq_)
          ##  msg = string(sym, "  in the '", eqg , "' equation ")
          ##  println(sym)
          ##  push!(errors, LinterWarning(errvalue, errtype, msg, loc, filename))
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
