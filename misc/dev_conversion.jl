using Dolo
using Dolo: transition
function hypeof(model)
    typ = typeof(model)
    otyp = eltype(model)
    ntyp = Float32
    s = replace(string(typ),string(otyp)=>"Float32")
    gtyp = ( Meta.parse(s) )
    return Union{eval(gtyp), typ}
end

# hypeof(model)
root_dir = pkgdir(Dolo)

# macro mdef(funtree)
#     println(funtree)
#     println(dump(funtree))
#     println("\nDONE\n")
#     return funtree
# end

model = include("$(root_dir)/examples/ymodels/rbc_mc.jl")


# macro mdef(funtree)
#     return funtree
# end

# @mdef function fun(a)
#     a = (;x=1,y=2)
#     (;x, y) = a
#     return (x*y)
# end



model32 = Dolo.convert_precision(Float32, model)



s_ = Dolo.calibrated(model32, :states)
x = Dolo.calibrated(model32, :controls)
m_ = Dolo.calibrated(model32, :exogenous)

s = [m_;s_]

dm32 = Dolo.discretize(model32)
Dolo.time_iteration(dm32, verbose=true)

Dolo.time_iteration(model32)

Dolo.transition(model32, s,x,m_)



function Dolo.transition(model::typeof(model), s::NamedTuple, x::NamedTuple, M::NamedTuple)
    
    (;δ, ρ) = model.calibration
    
    # Z = e.Z
    K = s.k * (1-δ) + x.i

    (;k=K,)  ## This is only the endogenous state

end


model_32 = Dolo.convert_precision(Float32, model)
# inherit all methods
transition(m::typeof(model_32), args...) = transition(model, args...)
arbitrage(m::typeof(model_32), args...) = arbitrage(model, args...)

Dolo.time_iteration(model_32)


dm32 = Dolo.discretize(model_32, Dict(:endo=>[1000]) )

typeof(dm32)

using Adapt
import oneAPI: oneArray

gpuArray = oneArray
interp_mode = :linear

wk0 = Dolo.time_iteration_workspace(dm32; interp_mode=interp_mode);
wk1 = Dolo.time_iteration_workspace(dm32; interp_mode=interp_mode);
wk_gpu = Dolo.time_iteration_workspace(dm32, dest=gpuArray; interp_mode=interp_mode);

@time sol1 = time_iteration(dm32, wk0; engine=:nothing, tol_η=1e-5, verbose=true, improve_wait=0, improve=false);
@time sol3 = time_iteration(dm32, wk1; engine=:cpu, tol_η=1e-5, verbose=true, improve_wait=0, improve=false);
@time sol2 = time_iteration(dm32, wk_gpu; engine=:gpu, tol_η=1e-5, verbose=true, improve_wait=0, improve=false);
