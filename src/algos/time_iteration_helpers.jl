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
# function euler_residuals_ti!(res::MSM{Point{n_x}}, model, dprocess::AbstractDiscretizedProcess,s::ListOfPoints{d}, x::MSM{Point{n_x}}, p::SVector, dr) where n_x where d

#     N = length(s)

#     for i in 1:length(res)
#         m = node(Point, dprocess, i)
#         for (w, M, j) in get_integration_nodes(Point, dprocess,i)
#             # Update the states
#             # TODO: replace views here
#             S = Dolo.transition(model, m, s, x.views[i], M, p)
#             X = dr(i, j, S)
#             res.views[i][:] += w*Dolo.arbitrage(model, m, s, x.views[i], M, S, X, p)
#         end
#     end
#     return res
# end

# function euler_residuals_ti(model, dprocess::AbstractDiscretizedProcess,s::ListOfPoints{d}, x::Vector{Vector{Point{n_x}}}, p::SVector, dr) where n_x where d
#     # res = deepcopy(x)
#     # for i_m=1:length(res)
#     #     res[i_m][:] *= 0.0
#     # end
#     xx = MSM(x)
#     res = copy(xx)
#     reset!(res)
#     euler_residuals_ti!(res, model, dprocess, s, xx, p, dr)
#     return vecvec(res)
# end

struct Euler{ID, n_m, n_s, n_x, Gx, Ge}

    model::AbstractModel
    p:: SVector  ## SVector
    grid::ProductGrid{Gx, Ge}
    grid_endo # obsolete
    grid_exo # obsolete
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

        n_m = length(model.symbols[:exogenous])
        n_s = length(model.symbols[:states])
        n_x = length(model.symbols[:controls])
        

        xxx = [ [dr0(j, Dolo.node(grid_endo, i)) for i=1:N_s] for j=1:N_m]

        x0 = MSM(xxx)

        s0 = Dolo.nodes(grid_endo)

        dr = Dolo.CachedDecisionRule(dprocess, grid_endo, xxx)
        cdr = Dolo.CachedDecisionRule(dprocess, grid_endo, xxx)


        m0 = SVector(model.calibration[:exogenous]...)

        grid = ProductGrid(grid_exo, grid_endo)
        Gx = typeof(grid.exo)
        Ge = typeof(grid.endo)


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

        ID = id(model)

        new{ID, n_m, n_s, n_x, Gx, Ge}(model, p, grid, grid_endo, grid_exo, dprocess, s0, x0, bounds, dr, cdr)

    end

end

repsvec(v1::Point{d1}, v2::Point{d}) where d1 where d= SVector(v1..., v2[d1+1:end]...)

function euler_residuals(F::Euler{id, n_m, n_s, n_x}, s::Vector{Point{n_s}}, x::MSM{Point{n_x}}, dr, 
    dprocess, parms::SVector; exo=nothing, keep_J_S=false, out=nothing) where id where n_s where n_x where n_m #, jres=nothing, S_ij=nothing)
    #
    # if set_dr ==true
    #   set_values!(dr,x)
    # end

    model = F.model

    if out === nothing
        res = zeros_like(x)::MSM{Point{n_x}}
    else
        res = out
    end

    res.data .*= 0.0 # just to be on the safe side...
    
    N_s = length(s) # Number of gris points for endo_var

    n_ms = length(x.views)  # number of exo states today
    n_mst = n_inodes(dprocess,1)  # number of exo states tomorrow
    d = length(s[1])

    # TODO: allocate properly...
    
    if keep_J_S
        jres = zeros((n_ms,n_mst,N_s,n_x,n_x))
        fut_S = zeros((n_ms,n_mst,N_s,n_s))
        J_ij = Vector{SMatrix{n_x,n_x,Float64,n_x*n_x}}[to_LOJ(jres[i,j,:,:,:]) for i=1:size(jres,1), j=1:size(jres,2)]
        S_ij = Vector{Point{d}}[to_LOP(fut_S[i,j,:,:]) for i=1:size(fut_S,1), j=1:size(fut_S,2)]
    end

    for i_ms in 1:n_ms
        m = node(Point{n_m}, dprocess,i_ms)
        if !(exo === nothing)
            m = repsvec(exo[1], m)   # z0
        end
        xx = x.views[i_ms]
        for I_ms in 1:n_mst
            M = inode(Point{n_m}, dprocess, i_ms, I_ms)
            if !(exo === nothing)
                M = repsvec(exo[2], M)   # z1
            end
            w = iweight(dprocess, i_ms, I_ms)
            S = transition(model, m, s, xx, M, parms)
            X = dr(i_ms, I_ms, S)
            if keep_J_S==true
                rr, rr_XM = arbitrage(model,Val{(0,6)},m,s,xx,M,S,X,parms)
                J_ij[i_ms,I_ms][:] = w*rr_XM
                S_ij[i_ms,I_ms][:] = S
                res.views[i_ms][:] += w*rr
            else
                rr = arbitrage(model, m, s, xx, M, S, X, parms)
                res.views[i_ms][:] .+= w*rr
            end
        end
    end

    # TODO: make this function type-stable
    if keep_J_S
        return res,J_ij,S_ij
    else
        return res
    end
end

function (F::Euler)(x0::MSM, x1::MSM, z0::Point{n_z}, z1::Point{n_z}; exo=nothing, kwargs...) where n_z
    F(x0, x1; exo=(z0,z1), kwargs...)
end

function (F::Euler)(x0::MSM, x1::MSM; exo=nothing, set_future=true, ignore_constraints=false, out=nothing, p=F.p) 

    if out === nothing
        rr = copy(x0)
    else
        rr = out
    end

    if set_future
        set_values!(F.dr, x1)
    end
    
    euler_residuals(F,F.s0,x0,F.dr,F.dprocess,p;exo=exo,keep_J_S=false,out=rr)
    
    if (F.bounds!==nothing) & !ignore_constraints
        lb, ub = F.bounds
        for n=1:length(rr.data)
            rr.data[n] = PhiPhi0(rr.data[n], x0.data[n], lb.data[n], ub.data[n])
        end
    end

    return rr

end

function df_A(F, z0, z1; set_future=false, inplace=false, ret_res=false, exo=nothing)

    if set_future
        set_values!(F.dr, z1)
    end

    if !inplace
        fun  = z->F(z, z1; set_future=false, ignore_constraints=true, exo=exo)
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
            f0.data[n] = z
            J.data[n]= z_f*J.data[n] + z_x
        end
    end

    if ret_res
        return (f0,J)
    else
        return J
    end
end


function df_B(F, x0, x1; set_future=false, exo=nothing)

    ddr_filt = F.cdr

    # res = deepcopy(z0)

    if set_future
        set_values!(F.dr, x1)
    end

    rr,J_ij,S_ij  = euler_residuals(F, F.s0, x0 , F.dr, F.dprocess, F.p; keep_J_S=true, exo=exo)

    if (F.bounds!==nothing)
        lb, ub = F.bounds
        PhiPhi!(rr.views, x1.views, lb.views, ub.views, J_ij)
    end

    L = LinearThing(J_ij, S_ij, ddr_filt)

    return L

end


function DiffFun2(fun, e::SVector{n_p}, ε=1e-6) where n_p

    # TODO: this could be done without allocation
    r0 = fun(e)
    n_x = length(r0.data[1])

    diffs = []
    for i=1:n_p
        ei = SVector([(j==i) ? e[i]+ε : e[j] for j=1:n_p ]...)
        d_i = (fun(ei)-r0)/ε
        push!(diffs, d_i)
    end

    res = [
        SMatrix{n_x, n_p, Float64, n_x*n_p}(cat([diffs[i].data[n] for i=1:n_p]...;dims=2))
        for n=1:length(r0.data)
    ]

    return MSM(res, r0.sizes)

end


function df_e(F, x0, x1, e1, e2; set_future=false)

    ε = 1e-8

    # diff w.r.t. e1
    # diff w.r.t. e2

    tt = F(x0,x1; exo=(e1,e2), set_future=false)

    # Awful hack
    fun1 = e -> F(x0,x1; exo=(e,e2))

    J1 = Dolo.DiffFun2(fun1, e1, 1e-8)
    
    fun2 = e -> F(x0,x1; exo=(e1,e))
    J2 = Dolo.DiffFun2(fun2, e2, 1e-8)

    return J1, J2

    
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