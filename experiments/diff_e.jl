using Dolo
import Dolo: Euler
model = Model("examples/models/rbc.yaml")

F = Euler(model);

function L_p(F, x0, x1, dp; ε=1e-8)

    p0 = F.p

    dp = SVector(dp...)

    u = maximum(abs, dp)


    dp_ε = dp/u*1e-8

    f0 = F(x0, x1)

    F.p = p0 + dp_ε

    f0_ϵ = F(x0, x1)

    F.p = p0

    f = (f0_ε-f0)/ε*u

    return f

end

