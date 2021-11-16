import Dolo: MSM

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

function (F::Euler)(x0::MSM, x1::MSM, set_future=true) 

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


function df_A(F, z0, z1; set_future=false)

    fun  = z->F(z, z1, false)

    J = Dolo.DiffFun(fun, z0)[2]

    return (J)

end


function df_B(F, z0, z1; set_future=false)

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

