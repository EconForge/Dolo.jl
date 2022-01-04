
import Base: -, \, +, /, *
import Base: copy
import Dolo: maxabs



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


norm(a::MSM) = maximum( u-> maximum(abs,u), a.data )
maxabs(a::MSM) = maximum( u-> maximum(abs,u), a.data )


function *(x::MSM, y::MSM)
    data = x.data .* y.data
    return MSM(data, x.sizes)
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


