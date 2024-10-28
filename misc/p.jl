using Dolo

model = include("../examples/ymodels/rbc_ar1.jl")

Dolo.discretize(model)


using StaticArrays

sg = Dolo.discretize(cs, [NN,NN,NN])

coeffs = [(i+1)^2 for i=1:d]
coeffs_2 = [(i+1)^3 for i=1:d]

fl(x) = SVector( sum(x .* coeffs), sum(x.* coeffs_2) )
fc(x) = SVector( sum(x .* coeffs + x.^2 .* coeffs), sum(x.* coeffs_2  + x.^3 .* coeffs) )
function testfun(dfun, sg)
    if sum(sum(dfun(e) for e in sg)) <0.0
        println("Oups")
    else
        nothing
    end
end


values = [fc(e) for e in sg]
gvec = Dolo.GVector(sg, values)


#linear 
values = [fl(e) for e in sg]
gvec = Dolo.GVector(sg, values)

dfun = Dolo.DFun(cs, gvec);

_values = [dfun(e) for e in sg]

d = _values - values
@test isapprox(values, _values)
testfun(dfun, sg)
@test 0== @allocated testfun(dfun, sg)


#cubic
using Test

dfun = Dolo.DFun(cs, gvec; interp_mode=:cubic)

_values = [dfun(e) for e in sg]

diff = _values - values

@test isapprox(values, _values)

testfun(dfun, sg)
@test 0== @allocated testfun(dfun, sg)



import Dolo.splines: prefilter_1!, prefilter_2!


v = rand(20)
w1 = copy(v)
w2 = copy(v)

prefilter_1!(w1)
prefilter_2!(w2)