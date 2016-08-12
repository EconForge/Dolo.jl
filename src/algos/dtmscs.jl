#=
+----------------------------------------------+
|                                              |
|             CompEcon algorithm               |
|                                              |
+----------------------------------------------+
=#
using CompEcon

type CEParams
    Ns
    Nm
    Nsf
    Nmf
    Ns_m
    Ns_mf
    Phi
    Phif
    Emat
    Qe
    Phi_m
    Phi_s
    Phi_mf
    sgrid
    sgri0
    sgridf
    mgrid
    mgrid0
    P
    S
    Sf
    spline_order
end

function CEParams(m::DTMSCCModel, spline_order, Nf)
end
