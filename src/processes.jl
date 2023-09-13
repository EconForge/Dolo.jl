### from QuantEcon
function gridmake!(out, arrays::Union{AbstractVector,AbstractMatrix}...)
    lens = Int[size(e, 1) for e in arrays]

    n = sum(_i -> size(_i, 2), arrays)
    l = prod(lens)
    @assert size(out) == (l, n)

    reverse!(lens)
    repititions = cumprod(vcat(1, lens[1:end-1]))
    reverse!(repititions)
    reverse!(lens)  # put lens back in correct order

    col_base = 0

    for i in 1:length(arrays)
        arr = arrays[i]
        ncol = size(arr, 2)
        outer = repititions[i]
        inner = floor(Int, l / (outer * lens[i]))
        for col_plus in 1:ncol
            row = 0
            for _1 in 1:outer, ix in 1:lens[i], _2 in 1:inner
                out[row+=1, col_base+col_plus] = arr[ix, col_plus]
            end
        end
        col_base += ncol
    end
    return out
end

@generated function gridmake(arrays::AbstractArray...)
    T = reduce(promote_type, eltype(a) for a in arrays)
    quote
        l = 1
        n = 0
        for arr in arrays
            l *= size(arr, 1)
            n += size(arr, 2)
        end
        out = Matrix{$T}(undef, l, n)
        gridmake!(out, arrays...)
        out
    end
end

function gridmake(t::Tuple)
    all(map(x -> isa(x, Integer), t)) ||
        error("gridmake(::Tuple) only valid when all elements are integers")
    gridmake(map(x->1:x, t)...)::Matrix{Int}
end

ckron(A::AbstractArray, B::AbstractArray) = kron(A, B)
ckron(arrays::AbstractArray...) = reduce(kron, arrays)

function rouwenhorst(N::Integer, ρ::Real, σ::Real, μ::Real=0.0)
    σ_y = σ / sqrt(1-ρ^2)
    p  = (1+ρ)/2
    ψ = sqrt(N-1) * σ_y
    m = μ / (1 - ρ)

    state_values, p = _rouwenhorst(p, p, m, ψ, N)
    (;P=p, V=state_values)
end

function _rouwenhorst(p::Real, q::Real, m::Real, Δ::Real, n::Integer)
    if n == 2
        return [m-Δ, m+Δ],  [p 1-p; 1-q q]
    else
        _, θ_nm1 = _rouwenhorst(p, q, m, Δ, n-1)
        θN = p    *[θ_nm1 zeros(n-1, 1); zeros(1, n)] +
             (1-p)*[zeros(n-1, 1) θ_nm1; zeros(1, n)] +
             q    *[zeros(1, n); zeros(n-1, 1) θ_nm1] +
             (1-q)*[zeros(1, n); θ_nm1 zeros(n-1, 1)]

        θN[2:end-1, :] ./= 2

        return range(m-Δ, stop=m+Δ, length=n), θN
    end
end

using Kronecker: kron

###


#### From Quantecon / miranda & Fackler

function qnwnorm(n::Int)
    maxit = 100
    pim4 = 1 / pi^(0.25)
    m = floor(Int, (n + 1) / 2)
    nodes = zeros(n)
    weights = zeros(n)

    z = sqrt(2n + 1) - 1.85575 * ((2n + 1).^(-1 / 6))

    for i = 1:m
        # Reasonable starting values for root finding
        if i == 1
            z = sqrt(2n + 1) - 1.85575 * ((2n + 1).^(-1 / 6))
        elseif i == 2
            z = z - 1.14 * (n.^0.426) ./ z
        elseif i == 3
            z = 1.86z + 0.86nodes[1]
        elseif i == 4
            z = 1.91z + 0.91nodes[2]
        else
            z = 2z + nodes[i - 2]
        end

        # root finding iterations
        it = 0
        pp = 0.0  # initialize pp so it is available outside while
        while it < maxit
            it += 1
            p1 = pim4
            p2 = 0.0

            for j = 1:n
                p3 = p2
                p2 = p1
                p1 = z .* sqrt(2 / j) .* p2 - sqrt((j - 1) / j) .* p3
            end

            # p1 now contains degree n Hermite polynomial
            # pp is derivative of p1 at the n'th zero of p1
            pp = sqrt(2n) .* p2
            z1 = z
            z = z1 - p1 ./ pp  # newton step

            if abs(z - z1) < 1e-14
                break
            end
        end

        if it >= maxit
            error("Failed to converge in qnwnorm")
        end

        nodes[n + 1 - i] = z
        nodes[i] = -z
        weights[i] = 2 ./ (pp .* pp)
        weights[n + 1 - i] = weights[i]
    end

    weights = weights ./ sqrt(pi)
    nodes = sqrt(2) .* nodes

    return nodes, weights
end

using LinearAlgebra: cholesky

function qnwnorm(n::Vector{Int}, mu::Vector, sig2::Matrix = Matrix(I, length(n), length(n)))
    n_n, n_mu = length(n), length(mu)

    if !(n_n == n_mu)
        error("n and mu must have same number of elements")
    end

    _nodes = Array{Vector{Float64}}(undef, n_n)
    _weights = Array{Vector{Float64}}(undef, n_n)

    for i in 1:n_n
        _nodes[i], _weights[i] = qnwnorm(n[i])
    end

    weights = ckron(_weights[end:-1:1]...)
    nodes = gridmake(_nodes...)::Matrix{Float64}

    new_sig2 = cholesky(sig2).U

    mul!(nodes, nodes, new_sig2)
    broadcast!(+, nodes, nodes, mu')

    return nodes, weights
end

struct MvNormal{names,n,n2}
    μ::SVector{n,Float64}
    Σ::SMatrix{n,n,Float64,n2}
end

MvNormal(names, μ::SVector{n,Float64}, Σ::SMatrix{n,n, Float64,n2}) where n where n2 = MvNormal{names, n, n2}(μ,Σ)
variables(mv::MvNormal{n}) where n = n

MvNormal(names, Σ::SMatrix{n,n, Float64,n2}) where n where n2 = MvNormal{names,n,n2}(zero(SVector{n,Float64}),Σ)

function MvNormal(names, μ::Vector, Σ::Matrix)
    p = size(Σ,1)
    v = SVector(μ...)
    M = SMatrix{p,p,Float64,p*p}(Σ)
    MvNormal(names, v, M)
end

function MvNormal(names, Σ::Matrix)
    MvNormal(names, zeros(size(Σ, 1)), Σ)
end


import Base: rand
import Distributions


function Base.rand(mv::MvNormal{names, n}) where names where n

    dis = Distributions.MvNormal(Matrix(mv.Σ))
    m = Distributions.rand(dis)
    SVector(m...)
end

function Base.rand(mv::MvNormal{names, 1}) where names
    return SVector(randn()*sqrt(mv.Σ[1,1]))
end

function discretize(mv::MvNormal, n::Int=5)
    nn = [n for i=1:size(mv.Σ,1)]
    μ = zeros(size(mv.Σ,1))
    x, w = qnwnorm(nn, μ, Matrix(mv.Σ))
    
    xm = SVector((SVector(x[i,:]...) for i=1:n)...)
    (;x=xm,w=SVector(w...))
end

discretize(mv::MvNormal, d::Dict) = length(d)>=1 ? discretize(mv,d[:n]) : discretize(mv)


struct VAR1{names,V,B}
    ρ::Float64
    Σ::B
end

VAR1(names, ρ, Σ) = VAR1{names, typeof(ρ), typeof(Σ)}(ρ, Σ)


function rand(var::VAR1{N,V,B}, m::SVector{d,Float64}) where N where V where B<:SMatrix{1,1,Float64,1} where d
    
    return SVector( m[1]*var.ρ + randn()*sqrt(var.Σ[1,1]) )

end

function rand(var::VAR1, m0::SVector{d,Float64}) where d
    ρ = var.ρ
    dis = Distributions.MvNormal(Matrix(var.Σ))
    m = ρ*m0 + rand(dis)
    SVector(m...)
end


function rand(var::VAR1, m0::NamedTuple)
    v = rand(var, SVector(m0...))
    return NamedTuple{variables(var)}(v)
end

variables(::VAR1{names, ρ, Σ}) where names where ρ where Σ = names

using Kronecker

function discretize(var::VAR1, n::Int=3)

    names = variables(var)

    # it would be good to have a special type of VAR1 process
    # which has a scalar autoregressive term
    ρ = var.ρ
    d = size(var.Σ, 1)
    # d = size(ρ, 1)
    # ρ = ρ[1, 1]
    # @assert maximum(abs, ρ.-Matrix(ρ*I,d,d))<1e-16


    if size(var.Σ, 1) == 1
        (;P, V) = rouwenhorst(n, ρ, sqrt(var.Σ[1,1]))
        mat = Matrix( Matrix(V')' )
        return MarkovChain(names, P,mat)
    end


    C = cholesky(var.Σ) # σ = L'*L

    if n isa Int
        NN = fill(n, d) # default: same number of nodes in each dimension
    else
        NN = n
    end

    components = [rouwenhorst(N, ρ, 1.0) for N in NN]
    # mc_components = [MarkovChain(names, mc.P, Matrix(mc.V')) for mc in components]
    # mc_prod = MarkovProduct(mc_components...)

    P = kronecker( (mc.P for mc in components)... )
    # V = cat( (mc.V' for mc in components)...; dims=1)

    # return V, C, components
    # return P, V
    n = size(P,1)
    V = SVector( (C.L * SVector(x...) for x in Iterators.product( (e.V for e in components)... ))... )
    return  MarkovChain(names, SMatrix{n,n}(P), V)

end

discretize(var::VAR1, d::Dict) = length(d)>=1 ? discretize(var, d[:n]) : discretize(var)

struct MarkovChain{names, d, d2, k}
    P::SMatrix{d,d,Float64,d2}
    Q::SVector{d, SVector{k, Float64}}
end

variables(mc::MarkovChain{names}) where names = names
ndims(mc::MarkovChain{names}) where names = length(names)

# MarkovChain(names, P, Q) = MarkovChain{names, typeof(P), typeof(Q)}(P,Q)

function MarkovChain(names, P::Matrix, Q::Matrix) 
    d = size(P,1)
    sm = SMatrix{d,d,Float64,d*d}(P)
    sv = SVector( tuple( ( SVector(Q[i,:]...) for i in 1:size(Q,1))...)  )
    MarkovChain{names,d,d^2,length(sv[1])}(sm,sv)
end

MarkovChain(names, P::SMatrix, Q::SVector{d,SVector{k,Float64}}) where d where k = MarkovChain{names, size(P,1), length(P), length(Q[1])}(P, Q) # TODO: specify type arguments
function MarkovChain(P::SMatrix, Q::SVector{d,SVector{k,Float64}}) where d where k
    names = tuple((Symbol(string("e", i)) for i=1:k)...)
    MarkovChain(names, P, Q)
end
   

MarkovProduct(mc::MarkovChain) = mc


discretize(mc::MarkovChain, args...; kwargs...) = mc

