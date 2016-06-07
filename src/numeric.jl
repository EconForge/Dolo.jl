# ------------------- #
# Numeric model types #
# ------------------- #
immutable DTCSCCModel{ID} <: ANM{ID,:dtcscc}
    symbolic::ASM
    calibration::ModelCalibration
    options::Dict{Symbol,Any}
    distribution::Dict{Symbol,Any}
    model_type::Symbol
    name::UTF8String
    filename::UTF8String
end

immutable DTMSCCModel{ID} <: ANM{ID,:dtmscc}
    symbolic::ASM
    calibration::ModelCalibration
    options::Dict{Symbol,Any}
    distribution::Dict{Symbol,Any}
    model_type::Symbol
    name::UTF8String
    filename::UTF8String
end

_numeric_mod_type{ID}(::ASM{ID,:dtcscc}) = DTCSCCModel{ID}
_numeric_mod_type{ID}(::ASM{ID,:dtmscc}) = DTMSCCModel{ID}

for TM in (:DTCSCCModel, :DTMSCCModel)
    @eval begin
        # # function type constructor
        # function $(TF)(sm::SymbolicModel; print_code::Bool=false)
        #     if model_type(sm) != model_type($TF)
        #         msg = string("Symbolic model is of type $(model_type(sm)) ",
        #                      "cannot create functions of type $($TF)")
        #         error(msg)
        #     end
        #     $(TF)([let
        #                eval(Dolo, compile_equation(sm, fld; print_code=print_code))
        #            end
        #            for fld in fieldnames($(TF))]...)
        # end

        Base.convert(::Type{SymbolicModel}, m::$(TM)) = m.symbolic

        # model type constructor
        function ($TM){ID}(sm::SymbolicModel{ID}; print_code::Bool=false)
            # compile all equations
            for eqn in keys(sm.equations)
                eval(Dolo, compile_equation(sm, eqn; print_code=print_code))
            end
            calib = ModelCalibration(sm)
            options = _to_Float64(eval_with(calib, deepcopy(sm.options)))

            # hack to parse normal transition matrix into a matrix instead of
            # a vector of vectors
            dist = eval_with(calib, deepcopy(sm.distribution))
            if haskey(dist, :Normal)
                n = length(calib[:shocks])
                dist[:Normal] = reshape(vcat(dist[:Normal]...), n, n)
            end

            dist = _to_Float64(dist)
            $(TM){ID}(sm, calib, options, dist, sm.model_type,
                      sm.name, sm.filename)
        end
    end
end

# ------------- #
# Other methods #
# ------------- #

filename(m::AbstractModel) = m.filename
name(m::AbstractModel) = m.name
