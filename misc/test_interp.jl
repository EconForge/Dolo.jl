using Test
using StaticArrays
using Dolo
import Dolo.splines: prefilter!, eval_UC_spline

d = 1


NN = 20


v = range(0,10;length=11)
v0 = [0.0; v ; 0.0]


prefilter!(v0)

vv = [SVector(e) for e in v]

res  = eval_UC_spline( (0.0, ), (1.0, ), (11,), v0, vv)







vars =  tuple( (Symbol("d$i") for i in 1:d)... )
pairs = (vars[i]=>[0.0, 1.0*i] for i in 1:d)
dvars = Dict( pairs )
cs = Dolo.CSpace(;
    pairs...
)

sg = Dolo.discretize(cs)

coeffs = [(i+1)^2 for i=1:d]
coeffs_2 = [(i+1)^3 for i=1:d]

fl(x) = SVector( sum(x .* coeffs), sum(x.* coeffs_2) )
fc(x) = SVector( sum(x .* coeffs + x.^2 .* coeffs), sum(x.* coeffs_2  + x.^3 .* coeffs) )
