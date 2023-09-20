using Dolo

model = include("examples/ymodels/rbc_iid.jl")


function get_defs(model)

    # could say where functions are defined (YAML, file, default)

    tm = typeof(model)
    defs = Dict()
    for fun in [Dolo.transition, Dolo.arbitrage]
        res = methods(fun)
        output = [r for r in res.ms if r.sig.types[2]==tm]
        name = output[1].name
        if length(output) == 1
            m = output[1]
            defs[name] = true
        else
            defs[name] = false
        end
    end

    return defs

end