__precompile__()

module Bruteforce_module

import Dolo
using AxisArrays

export euler_residuals, SerialDifferentiableFunction, swaplines, divide, substract,
       invert, ssmul, destack0, d_filt_dx, invert_jac

export ITI_function
# export ImprovedTimeIterationResult, Base.show


include("bruteforce_help.jl")
include("ITI_function.jl")

end
