using Dolo
using Test

d = 2

NN = 20

vars =  tuple( (Symbol("d$i") for i in 1:d)... )
sizes = tuple( (NN+i+1 for i in 1:d)...)
args = tuple( ( (0.0, 1.0*i, sizes[i]) for i in 1:d)... )
cg = Dolo.CGrid( args )

@test isbits(cg)

@test Dolo.size(cg) == sizes


# check that iterator for is non-allocating
f(cg,v) = begin e = sum( sum( e for e in cg ) ); v[1] = e ; nothing end
res = zeros(1)
f(cg, res)
@test 0==(@allocated f(cg, res))

# Dolo.enum
# returns QP object
# ...

# check that enum for is non-allocating
function fun(cg::Dolo.CGrid{d}) where d
    if sum( sum( (sum(loc)*val) for (;loc, val) in Dolo.enum(cg) ) ) < -Inf
        print("")
    end
end

fun(cg)
@test 0==(@allocated fun(cg))
@time fun(cg)


println((@allocated fun(cg)))
