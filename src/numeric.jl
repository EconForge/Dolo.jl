# ------------------- #
# Numeric model types #
# ------------------- #

immutable Options{TD<:AbstractDistribution, TG<:AbstractGrid}
    grid::TG
    distribution::TD
    other::Dict{Symbol,Any}  # TODO: shouldn't need. Just keeps stuff around
end

Options(;grid=nothing, distribution=nothing, other=Dict{Symbol,Any}()) =
    Options(grid, distribution, other)

function Options(sm::AbstractSymbolicModel, calib::ModelCalibration)
    # numericize options
    _options = eval_with(calib, deepcopy(sm.options))

    _opts = Dict{Symbol,Any}()
    other = Dict{Symbol,Any}()

    # now construct Options object
    for k in keys(_options)
        data = pop!(_options, k)
        if k == :grid
            _opts[:grid] = _build_grid(data)
        elseif k == :distribution
            _opts[:distribution] = _build_dist(data, calib)
        else
            other[k] = data
        end
    end

    Options(;_opts..., other=other)
end

# TODO: given that fields are exactly the same should we make just a single
#       NumericModel and distinguish between DTCSCC and DTMSCC via type params?
immutable NumericModel{ID,kind} <: ANM{ID,kind}
    symbolic::SymbolicModel{ID,kind}
    calibration::ModelCalibration
    options::Options
    model_type::Symbol
    name::UTF8String
    filename::UTF8String
end

typealias DTCSCCModel{ID} NumericModel{ID,:dtcscc}
typealias DTMSCCModel{ID} NumericModel{ID,:dtmscc}
typealias DynareModel{ID} NumericModel{ID,:dynare}

_numeric_mod_type{ID}(::ASM{ID,:dtcscc}) = DTCSCCModel{ID}
_numeric_mod_type{ID}(::ASM{ID,:dtmscc}) = DTMSCCModel{ID}
_numeric_mod_type{ID}(::ASM{ID,:dynare}) = DynareModel{ID}

function compile_equation(sm::SymbolicModel, func_nm::Symbol; print_code::Bool=false)
    # extract spec from recipe
    spec = RECIPES[model_type(sm)][:specs][func_nm]

    # get expressions from symbolic model
    exprs = sm.equations[func_nm]

    numeric_mod = _numeric_mod_type(sm)

    bang_func_nm = Symbol(string(func_nm), "!")

    if length(exprs) == 0
        msg = "Model did not specify functions of type $(func_nm)"
        code = quote
            function $(func_nm)(::$(numeric_mod), args...)
                error($msg)
            end

            function $(bang_func_nm)(::$(numeric_mod), args...)
                error($msg)
            end
        end
    else
        code = make_method(FunctionFactory(sm, func_nm))
    end
    print_code && println(code)
    code
end

Base.convert(::Type{SymbolicModel}, m::NumericModel) = m.symbolic

function NumericModel{ID,kind}(sm::SymbolicModel{ID,kind}; print_code::Bool=false)
    # compile all equations
    for eqn in keys(sm.equations)
        eval(Dolo, compile_equation(sm, eqn; print_code=print_code))
    end

    # get numerical calibration and options
    calib = ModelCalibration(sm)
    options = Options(sm, calib)

    NumericModel(sm, calib, options, sm.model_type, sm.name, sm.filename)
end

# ------------- #
# Other methods #
# ------------- #

filename(m::AbstractModel) = m.filename
name(m::AbstractModel) = m.name
