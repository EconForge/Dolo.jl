d = 1

NN = 20

vars =  tuple( (Symbol("d$i") for i in 1:d)... )
sizes = tuple( (NN+i+1 for i in 1:d)...)
args = tuple( ( (0.0, 1.0*i, sizes[i]) for i in 1:d)... )
cg = Dolo.CGrid( args )

@test isbits(cg)

@test Dolo.size(cg) == sizes


# check that iterator for is non-allocating


# check that 'iterate' is non-allocating
function fun1(cg::Dolo.CGrid{d}) where d
    if sum( sum( e for e in cg) ) < -Inf
        print("This should never be ran.")
    end
end
fun1(cg)
@test 0==(@allocated fun1(cg))

# Dolo.enum
# returns QP object
# ...

# check that 'enum' for is non-allocating
function fun(cg::Dolo.CGrid{d}) where d
    if sum( sum( sum(q.loc...)*(q.val) for q in Dolo.enum(cg) ) ) < 0.0
        print("")
    end
end

fun(cg)
@test 0==(@allocated fun(cg))
