### What is the right way to do that ?
using RecipesBase

# Cartesian x scalar values
@recipe function f(φ::DFun{A,B}, n::Integer = 100) where A where B<:GArray{G,V} where G<:CGrid where V<:AbstractVector{Float64}
    
    xlabel --> string(dims(φ.domain)[1])
    ylabel --> string(vars(φ))
    
    dom = φ.domain
    xvec = range(dom.min[1], dom.max[1], n)  
    yvec = φ.(xvec)

    @series begin
        xx = [e[1] for e in φ.values.grid]       
        yy = [e for e in φ.values]
        markershape := :circle
        linealpha := 0.0
        (xx,yy)
    end

    return (xvec, yvec)

end