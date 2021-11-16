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
    res = deepcopy(xx)
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
    bounds::Union{Nothing, Tuple{MSM{SVector{n_x, Float64}}}}
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
        
        N_m = Dolo.n_nodes(grid_exo)
        N_s = Dolo.n_nodes(grid_endo)

        n_s = length(model.symbols[:states])
        n_x = length(model.symbols[:controls])
        

        xxx = [ [dr0(j, Dolo.node(grid_endo, i)) for i=1:N_s] for j=1:N_m]

        x0 = MSM(xxx)

        s0 = Dolo.nodes(grid_endo)

        dr = Dolo.CachedDecisionRule(dprocess, grid_endo, xxx)
        cdr = Dolo.CachedDecisionRule(dprocess, grid_endo, xxx)


        Gx = typeof(grid_exo)
        Ge = typeof(grid_endo)

        if :arbitrage_lb in keys(model.factories)
            lb = [Dolo.arbitrage_lb(model, m, s, p) for m in Dolo.nodes(grid_exo)]
            ub = [Dolo.arbitrage_ub(model, m, s, p) for m in Dolo.nodes(grid_exo)]
            bounds = (MSM(lb), MSM(ub))
        else
            bounds = nothing
        end

        new{n_s, n_x, Gx, Ge}(model, p, grid_endo, grid_exo, dprocess, s0, x0, bounds, dr, cdr)

    end

end

function (F::Euler)(x0::MSM, x1::MSM, set_future=true, ignore_constraints=False) 

    res = deepcopy(x0)

    if set_future
        set_values!(F.dr, x1)
    end
    
    rr =   euler_residuals(F.model,F.s0,x0,F.dr,F.dprocess,F.p;keep_J_S=false)

    return rr

end

import Base: -, \, +, /, *

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


function df_A(F, z0, z1; set_future=false, ignore_constraints=False)

    fun  = z->F(z, z1, false)

    J = Dolo.DiffFun(fun, z0)[2]

    return (J)

end


function df_B(F, z0, z1; set_future=false, ignore_constraints=False)

    ddr_filt = F.cdr

    # res = deepcopy(z0)

    if set_future
        set_values!(F.dr, z1)
    end

    _,J_ij,S_ij  =   euler_residuals(F.model, F.s0, z0 , F.dr, F.dprocess, F.p; keep_J_S=true)

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

function mult!(L::LinearThing,x::Number)
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

