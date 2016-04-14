using Dolo

# 3-dimensional space
a = [0.1, -0.2, 0.0]
b = [0.5, 1.0, 10.0]
dims = [10,10,10]
grid = mlinspace(a,b,dims)
ngrid = mlinspace(a,b,[20,20,20])

# values to interpolate
values = rand(prod(dims),3)

dr = DecisionRule(a,b,dims)
set_values(dr, values)

# vectorized evaluation
evaluate(dr, ngrid)

# evaluate at one point
evaluate(dr, slice(ngrid,1,:) )



### Multivalued Decision rules
mvalues = rand(4,size(grid,1), 3)

mdr = MixedDecisionRule(4,a,b,dims)
set_values(mdr, mvalues)

# evaluate at one point
evaluate(mdr, 2, slice(ngrid,1,:))

# vectorized evaluation
evaluate(mdr, 2, ngrid)
