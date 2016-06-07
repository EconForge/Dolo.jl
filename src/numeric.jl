# ------------------- #
# Numeric model types #
# ------------------- #
# TODO: decide if we really want all 8 type params. It doesn't hurt, it just
#       looks funny
immutable DTCSCCfunctions{T1,T2,T3,T4,T5,T6,T7,T8}
    arbitrage::T1
    transition::T2
    auxiliary::T3
    value::T4
    expectation::T5
    direct_response::T6
    controls_lb::T7
    controls_ub::T8
end

# NOTE: a type parameter is needed so the `functions` field is not abstract
#       when instances of this type are actually created.
immutable DTCSCCModel{_T<:DTCSCCfunctions} <: ANM
    symbolic::SymbolicModel
    functions::_T
    calibration::ModelCalibration
    options::Dict{Symbol,Any}
    distribution::Dict{Symbol,Any}
    model_type::Symbol
    name::UTF8String
    filename::UTF8String
end

immutable DTMSCCfunctions{T1,T2,T3,T4,T5,T6,T7,T8,T9}
    arbitrage::T1
    transition::T2
    auxiliary::T3
    value::T4
    expectation::T5
    direct_response::T6
    controls_lb::T7
    controls_ub::T8
    arbitrage_2::T9
end

immutable DTMSCCModel{_T<:DTMSCCfunctions} <: ANM
    symbolic::SymbolicModel
    functions::_T
    calibration::ModelCalibration
    options::Dict{Symbol,Any}
    distribution::Dict{Symbol,Any}
    model_type::Symbol
    name::UTF8String
    filename::UTF8String
end

for (TF, TM, ms) in [(:DTCSCCfunctions, :DTCSCCModel, :(:dtcscc)),
                     (:DTMSCCfunctions, :DTMSCCModel, :(:dtmscc))]
    @eval begin
        model_type(::$(TM)) = $ms
        model_type(::Type{$(TM)}) = $ms
        model_type(::$(TF)) = $ms
        model_type(::Type{$(TF)}) = $ms

        # function type constructor
        function $(TF)(sm::SymbolicModel; print_code::Bool=false)
            if model_type(sm) != model_type($TF)
                msg = string("Symbolic model is of type $(model_type(sm)) ",
                             "cannot create functions of type $($TF)")
                error(msg)
            end
            $(TF)([let
                       eval(Dolo, compile_equation(sm, fld; print_code=print_code))
                   end
                   for fld in fieldnames($(TF))]...)
        end
        Base.convert(::Type{SymbolicModel}, m::$(TM)) = m.symbolic

        # model type constructor
        function ($TM)(sm::SymbolicModel; print_code::Bool=false)
            if model_type(sm) != model_type($TM)
                msg = string("Symbolic model is of type $(model_type(sm)) ",
                             "cannot create model of type $($TM)")
                error(msg)
            end
            calib = ModelCalibration(sm)
            options = _to_Float64(eval_with(calib, deepcopy(sm.options)))
            dist = eval_with(calib, deepcopy(sm.distribution))
            # hack to parse normal transition matrix into a matrix instead of
            # a vector of vectors
            if haskey(dist, :Normal)
                n = length(calib[:shocks])
                dist[:Normal] = reshape(vcat(dist[:Normal]...), n, n)
            end
            dist = _to_Float64(dist)
            funcs = $(TF)(sm; print_code=print_code)
            $(TM)(sm, funcs, calib, options, dist, sm.model_type,
                  sm.name, sm.filename)
        end
    end
end

# ------------- #
# Other methods #
# ------------- #

model_type(m::AbstractModel) = m.model_type
filename(m::AbstractModel) = m.filename
name(m::AbstractModel) = m.name
