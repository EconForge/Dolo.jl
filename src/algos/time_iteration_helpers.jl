import Dolo: MSM



"""
Computes the residuals of the arbitrage equations. The general form of the arbitrage equation is

    `0 = E_t [f(m, s, x, M, S, X; p)]`

where `m` are current exogenous variables, `s` are current states,
`x` are current controls, `M` are next period's exogenous variables, `S` are next period's states, `X` are next period's controls, and `p` are the model parameters. This function evaluates the right hand side of the arbitrage equation for the given inputs.

If the list of current controls `x` is provided as a two-dimensional array (`ListOfPoints`), it is transformed to a one-dimensional array (`ListOfListOfPoints`).


# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `dprocess`: Discretized exogenous process.
* `s::ListOfPoints`: List of state variable values.
* `x::ListOfListOfPoints`: List of control variable values associated with each exogenous shock.
* `p::Vector{Float64}`: Model parameters.
* `dr`: Current guess for the decision rule.
# Returns
* `res`: Residuals of the arbitrage equation associated with each exogenous shock.
"""
function euler_residuals_ti!(res::MSM{Point{n_x}}, model, dprocess::AbstractDiscretizedProcess,s::ListOfPoints{d}, x::MSM{Point{n_x}}, p::SVector, dr) where n_x where d

    N = length(s)

    for i in 1:length(res)
        m = node(Point, dprocess, i)
        for (w, M, j) in get_integration_nodes(Point, dprocess,i)
            # Update the states
            # TODO: replace views here
            S = Dolo.transition(model, m, s, x.views[i], M, p)
            X = dr(i, j, S)
            res.views[i][:] += w*Dolo.arbitrage(model, m, s, x.views[i], M, S, X, p)
        end
    end
    return res
end

function euler_residuals_ti(model, dprocess::AbstractDiscretizedProcess,s::ListOfPoints{d}, x::Vector{Vector{Point{n_x}}}, p::SVector, dr) where n_x where d
    # res = deepcopy(x)
    # for i_m=1:length(res)
    #     res[i_m][:] *= 0.0
    # end
    xx = MSM(x)
    res = copy(xx)
    reset!(res)
    euler_residuals_ti!(res, model, dprocess, s, xx, p, dr)
    return vecvec(res)
end


struct Euler{n_s, n_x, Gx, Ge}

    model::AbstractModel
    p:: SVector  ## SVector
    grid_endo
    grid_exo
    dprocess
    s0::ListOfPoints{n_s}
    x0::MSM{SVector{n_x, Float64}}
    bounds::Union{Nothing,Tuple{MSM{SVector{n_x, Float64}},MSM{SVector{n_x, Float64}}}}
    dr::CachedDecisionRule
    cdr::CachedDecisionRule

    function Euler(model; 
        discretization=Dict(), 
        interpolation=:cubic,
        dr0=Dolo.ConstantDecisionRule(model),
        ignore_constraints=false
    )

        p = SVector(model.calibration[:parameters]...)
        
        grid, dprocess = discretize(model; discretization...)
        grid_endo = grid.endo
        grid_exo = grid.exo

        x = SVector(model.calibration[:controls]...)
        
        N_m = max(Dolo.n_nodes(grid_exo),1)
        N_s = Dolo.n_nodes(grid_endo)

        n_s = length(model.symbols[:states])
        n_x = length(model.symbols[:controls])
        

        xxx = [ [dr0(j, Dolo.node(grid_endo, i)) for i=1:N_s] for j=1:N_m]

        x0 = MSM(xxx)

        s0 = Dolo.nodes(grid_endo)

        dr = Dolo.CachedDecisionRule(dprocess, grid_endo, xxx)
        cdr = Dolo.CachedDecisionRule(dprocess, grid_endo, xxx)


        m0 = SVector(model.calibration[:exogenous]...)

        Gx = typeof(grid_exo)
        Ge = typeof(grid_endo)

        if (:controls_lb in keys(model.factories)) & !ignore_constraints
            if Dolo.n_nodes(grid_exo)>0
                lb = [Dolo.controls_lb(model, m, s0, p) for m in Dolo.nodes(grid_exo)]
                ub = [Dolo.controls_ub(model, m, s0, p) for m in Dolo.nodes(grid_exo)]
            else
                lb = [Dolo.controls_lb(model, m0, s0, p) ]
                ub = [Dolo.controls_ub(model, m0, s0, p) ]
            end
            bounds = (MSM(lb), MSM(ub))
        else
            bounds = nothing
        end

        new{n_s, n_x, Gx, Ge}(model, p, grid_endo, grid_exo, dprocess, s0, x0, bounds, dr, cdr)

    end

end

function (F::Euler)(x0::MSM, x1::MSM; set_future=true, ignore_constraints=false, out=nothing) 

    if out === nothing
        rr = copy(x0)
    else
        rr = out
    end

    if set_future
        set_values!(F.dr, x1)
    end
    
    # euler_residuals(F.model,F.s0,x0,F.dr,F.dprocess,F.p;keep_J_S=false, out=rr)
    euler_residuals(F.model,F.s0,x0,F.dr,F.dprocess,F.p;keep_J_S=false, out=rr)
    # euler_residuals_noalloc(F.model,F.s0,x0,F.dr,F.dprocess,F.p;keep_J_S=false, out=rr)

    if (F.bounds!==nothing) & !ignore_constraints
        lb, ub = F.bounds
        for n=1:length(rr.data)
            rr.data[n] = PhiPhi0(rr.data[n], x0.data[n], lb.data[n], ub.data[n])
        end
    end

    return rr

end

import Base: -, \, +, /, *
import Base: copy

function copy(m::MSM)
    dd = deepcopy(m.data)
    return MSM(dd, m.sizes)
end

function -(a::MSM)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = -a.data
    return MSM(data, sizes)
end


function -(a::MSM, b::MSM)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data .- b.data
    return MSM(data, sizes)
end



function +(a::MSM, b::MSM)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data .+ b.data
    return MSM(data, sizes)
end


function \(a::MSM, b::MSM)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data .\ b.data
    return MSM(data, sizes)
end

function /(a::MSM, b::MSM)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data ./ b.data
    return MSM(data, sizes)
end

function /(a::MSM, b::Number)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data ./ b
    return MSM(data, sizes)
end

function *(a::MSM, b::Number)
    N = length(a.data)
    sizes = [length(e) for e in a.views]
    # c = zeros_like(b)
    data = a.data .* b
    return MSM(data, sizes)
end



import Dolo: maxabs

norm(a::MSM) = maximum( u-> maximum(abs,u), a.data )
maxabs(a::MSM) = maximum( u-> maximum(abs,u), a.data )


function df_A(F, z0, z1; set_future=false, inplace=false)

    if set_future
        set_values!(F.dr, z1)
    end

    if !inplace
        fun  = z->F(z, z1; set_future=false, ignore_constraints=true)
        f0,J = Dolo.DiffFun(fun, z0, 1e-8)

    else
        n_x = length(z0.data[1])

        f0 = z0*0
        fi = z0*0
        xi = z0*0

        N = length(xi.data)
        cata = zeros(SMatrix{n_x, n_x, Float64, n_x*n_x}, N)
        J = MSM(cata, f0.sizes)
        out = f0, fi, xi, J

        fun!  = (r,z)->F(z, z1; set_future=false, ignore_constraints=true, out=r)

        f00, J0 = Dolo.DiffFun!(fun!, z0, 1e-8, out)

        @assert f00 === f0
        @assert J === J0

    end
    if (F.bounds!==nothing)
        lb, ub = F.bounds
        for n=1:length(f0.data)
            z,z_f,z_x  = PhiPhi(f0.data[n], z0.data[n], lb.data[n], ub.data[n])
            ### TODO: check PhiPhi and avoid -1 multiplication
            f0.data[n] = z
            J.data[n]= z_f*J.data[n] + z_x
        end
    end

    return (J)

end


function df_B(F, z0, z1; set_future=false)

    ddr_filt = F.cdr

    # res = deepcopy(z0)

    if set_future
        set_values!(F.dr, z1)
    end

    rr,J_ij,S_ij  =   euler_residuals(F.model, F.s0, z0 , F.dr, F.dprocess, F.p; keep_J_S=true)

    if (F.bounds!==nothing)
        lb, ub = F.bounds
        PhiPhi!(rr.views, z1.views, lb.views, ub.views, J_ij)
    end

    L = LinearThing(J_ij, S_ij, ddr_filt)

    return L

end


function *(x::MSM, y::MSM)
    data = x.data .* y.data
    return MSM(data, x.sizes)
end


function prediv!(L::LinearThing,x::MSM)
    for i=1:size(L.M_ij,1)
        N = length(x.views[i])
        for j=1:size(L.M_ij,2)
            for n=1:N
                L.M_ij[i,j][n] = x.views[i][n] \ L.M_ij[i,j][n]
            end
        end

    end
    return L
end

function premult!(L::LinearThing,x::MSM)
    for i=1:size(L.M_ij,1)
        N = length(x.views[i])
        for j=1:size(L.M_ij,2)
            for n=1:N
                L.M_ij[i,j][n] = x.views[i][n] * L.M_ij[i,j][n]
            end
        end

    end
    return L
end

function mult!(L::LinearThing,x::T) where T<:Number
    for i=1:size(L.M_ij,1)
        for j=1:size(L.M_ij,2)
            L.M_ij[i,j][:] *= x
        end

    end
end

function invert!(x::MSM)
    for i=1:length(x.data)
        x.data[i] = inv(x.data[i] )
    end
end



function ldiv(a::MSM{S}, b::MSM{T}) where S where T
    r = a.data .\ b.data
    MSM(r, a.sizes)
end

