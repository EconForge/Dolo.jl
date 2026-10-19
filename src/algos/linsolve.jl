import LinearAlgebra: dot
using TimerOutputs: TimerOutput, @timeit, disable_timer!

# default timer: disabled, so that @timeit has negligible overhead
const NO_TIMER = let t = TimerOutput(); disable_timer!(t); t end

# Linear solvers used by accelerated time iteration and Newton.
# Operators are passed as in-place functions A!(w, v) acting on flat vectors.

# Vectors of the size of the problem (r, w, aux, V) are stored where the data is (e.g. a GPU);
# the small quantities of the Krylov space (H, y, cs, sn, g) on the CPU.
struct GMRESWorkspace{Tf, TV<:AbstractVector{Tf}, TM<:AbstractMatrix{Tf}}
    r::TV
    w::TV
    aux::TV     # scratch vector for operators
    V::TM       # Krylov basis
    H::Matrix{Tf}   # Hessenberg matrix (triangularized in place)
    y::Vector{Tf}
    cs::Vector{Tf}  # Givens cosines
    sn::Vector{Tf}  # Givens sines
    g::Vector{Tf}   # rotated right-hand side
end

import Adapt
Adapt.adapt_structure(to, wk::GMRESWorkspace) = GMRESWorkspace(
    Adapt.adapt(to, wk.r), Adapt.adapt(to, wk.w), Adapt.adapt(to, wk.aux), Adapt.adapt(to, wk.V),
    wk.H, wk.y, wk.cs, wk.sn, wk.g
)

function gmres_workspace(Tf, n, m)
    GMRESWorkspace(
        zeros(Tf, n),
        zeros(Tf, n),
        zeros(Tf, n),
        zeros(Tf, n, m+1),
        zeros(Tf, m+1, m),
        zeros(Tf, m),
        zeros(Tf, m),
        zeros(Tf, m),
        zeros(Tf, m+1),
    )
end

krylov_dim(wk::GMRESWorkspace) = size(wk.H, 2)

"""
    gmres!(x, A!, b, wk; tol=1e-8, tmr=NO_TIMER)

Solve `A*x = b` with (unrestarted) GMRES, using the content of `x` as initial guess.
`A!(w, v)` computes `w = A*v` in place. The maximum number of iterations is the
dimension of the workspace `wk` (see `gmres_workspace`). The tolerance is on the
absolute residual. Returns `(niter, converged)`.
"""
function gmres!(x::AbstractVector{Tf}, A!, b::AbstractVector{Tf}, wk::GMRESWorkspace; tol=1e-8, tmr=NO_TIMER) where Tf

    (; r, w, V, H, y, cs, sn, g) = wk
    m = krylov_dim(wk)

    A!(w, x)
    @. r = b - w
    β = norm(r)
    if β < tol
        return 0, true
    end

    @views @. V[:,1] = r / β
    g .= zero(Tf)
    g[1] = β

    k = 0
    for j=1:m

        @timeit tmr "compute A.v" A!(w, view(V, :, j))

        # Arnoldi (modified Gram-Schmidt)
        for i=1:j
            H[i,j] = dot(view(V, :, i), w)
            @views @. w -= H[i,j] * V[:,i]
        end
        h = norm(w)
        H[j+1,j] = h

        # apply previous Givens rotations to the new column
        for i=1:j-1
            t = cs[i] * H[i,j] + sn[i] * H[i+1,j]
            H[i+1,j] = -sn[i] * H[i,j] + cs[i] * H[i+1,j]
            H[i,j] = t
        end

        # new rotation eliminating H[j+1,j]
        ρ = sqrt(H[j,j]^2 + H[j+1,j]^2)
        cs[j] = H[j,j] / ρ
        sn[j] = H[j+1,j] / ρ
        H[j,j] = ρ
        H[j+1,j] = zero(Tf)

        g[j+1] = -sn[j] * g[j]
        g[j] = cs[j] * g[j]

        k = j

        # lucky breakdown: the Krylov space is invariant and the solution is exact
        if abs(g[j+1]) < tol || h <= eps(Tf)*β
            break
        end

        if j < m
            @views @. V[:,j+1] = w / h
        end

    end

    # back substitution on the k×k triangular system
    for i=k:-1:1
        y[i] = g[i]
        for l=i+1:k
            y[i] -= H[i,l] * y[l]
        end
        y[i] /= H[i,i]
    end

    # y is on the CPU, V where the data is
    yk = similar(x, k)
    copyto!(yk, y[1:k])
    @views mul!(x, V[:, 1:k], yk, one(Tf), one(Tf))

    return k, abs(g[k+1]) < tol

end
