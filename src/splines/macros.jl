using StaticArrays

# TODO : for some reason, this doesn'work !
#
# const Ad_ = @SMatrix [
#     -1.0/6.0  3.0/6.0 -3.0/6.0 1.0/6.0;
#      3.0/6.0 -6.0/6.0  0.0/6.0 4.0/6.0;
#     -3.0/6.0  3.0/6.0  3.0/6.0 1.0/6.0;
#      1.0/6.0  0.0/6.0  0.0/6.0 0.0/6.0
#  ]
#  const dAd_ = @SMatrix [
#     0.0 -0.5  1.0 -0.5;
#     0.0  1.5 -2.0  0.0;
#     0.0 -1.5  1.0  0.5;
#     0.0  0.5  0.0  0.0
#  ]
#  const d2Ad_ = @SMatrix [
#     0.0 0.0 -1.0  1.0;
#     0.0 0.0  3.0 -2.0;
#     0.0 0.0 -3.0  1.0;
#     0.0 0.0  1.0  0.0
#  ]

const Ad64 = [
   [-1.0/6.0  3.0/6.0 -3.0/6.0 1.0/6.0];
   [ 3.0/6.0 -6.0/6.0  0.0/6.0 4.0/6.0];
   [-3.0/6.0  3.0/6.0  3.0/6.0 1.0/6.0];
   [ 1.0/6.0  0.0/6.0  0.0/6.0 0.0/6.0]
]

const dAd64 = [
   [ 0.0 -0.5  1.0 -0.5];
   [ 0.0  1.5 -2.0  0.0];
   [ 0.0 -1.5  1.0  0.5];
   [ 0.0  0.5  0.0  0.0]
]

const d2Ad64 = [
   [ 0.0 0.0 -1.0  1.0];
   [ 0.0 0.0  3.0 -2.0];
   [ 0.0 0.0 -3.0  1.0];
   [ 0.0 0.0  1.0  0.0]
]

const Ad32 = [
   [-1.0f0/6.0f0  3.0f0/6.0f0 -3.0f0/6.0f0 1.0f0/6.0f0];
   [ 3.0f0/6.0f0 -6.0f0/6.0f0  0.0f0/6.0f0 4.0f0/6.0f0];
   [-3.0f0/6.0f0  3.0f0/6.0f0  3.0f0/6.0f0 1.0f0/6.0f0];
   [ 1.0f0/6.0f0  0.0f0/6.0f0  0.0f0/6.0f0 0.0f0/6.0f0]
]

const dAd32 = [
   [ 0.0f0 -0.5f0  1.0f0 -0.5f0];
   [ 0.0f0  1.5f0 -2.0f0  0.0f0];
   [ 0.0f0 -1.5f0  1.0f0  0.5f0];
   [ 0.0f0  0.5f0  0.0f0  0.0f0]
]

const d2Ad32 = [
   [ 0.0f0 0.0f0 -1.0f0  1.0f0];
   [ 0.0f0 0.0f0  3.0f0 -2.0f0];
   [ 0.0f0 0.0f0 -3.0f0  1.0f0];
   [ 0.0f0 0.0f0  1.0f0  0.0f0]
]

U(s,i) = Symbol(string(s,i))
U(s,i,j) = Symbol(string(s,i,"_",j))

function create_Phi(d, extrap, diff; Tf=Float64)
    if Tf==Float64
        Ad = Ad64
        dAd = dAd64
    elseif Tf==Float32
        Ad = Ad32
        dAd = dAd32
    end
    lines = []
    for i=1:d
        block = []
        rhs_1 = U("tp", i,1)
        rhs_2 = U("tp", i,2)
        rhs_4 = U("tp", i,4)
        rhs_3 = U("tp", i,3)
        if extrap == "none"
            for j=1:4
                eq = :($(U("Phi_",i,j)) = ($(Ad[j,1])*$rhs_1 + $(Ad[j,2])*$rhs_2 + $(Ad[j,3])*$rhs_3 + $(Ad[j,4])*$rhs_4) )
                push!(lines,eq)
            end
        elseif extrap == "natural"
            block = quote
                if $(U("t",i))<0
                    $( [ :($(U("Phi_",i,j)) = ($(dAd[j,4])*$rhs_3 + $(Ad[j,4])) ) for j=1:4 ]...)
                elseif $(U("t",i))>1
                    $( [ :($(U("Phi_",i,j)) = (3*$(Ad[j,1]) + 2*$(Ad[j,2]) + $(Ad[j,3]))*($rhs_3-1) + ($(Ad[j,1])+$(Ad[j,2])+$(Ad[j,3])+$(Ad[j,4])) ) for j=1:4 ]...)
                else
                    $( [ :($(U("Phi_",i,j)) = ($(Ad[j,1])*$rhs_1 + $(Ad[j,2])*$rhs_2 + $(Ad[j,3])*$rhs_3 + $(Ad[j,4])*$rhs_4) ) for j=1:4 ]...)
                end
            end
            # for l in block.args
            push!(lines,block.args...)
        end
        if diff
            for j=1:4
                eq = :($(U("dPhi_",i,j)) = ($(dAd[j,1])*$rhs_1 + $(dAd[j,2])*$rhs_2 + $(dAd[j,3])*$rhs_3 + $(dAd[j,3])*$rhs_4 )*$(U("dinv",i)) )
                push!(lines,eq)
            end
        end
    end
    return lines
end

function tensor_prod(symbs, inds, add_index=false)
    if length(symbs)==0
        subscripts = [:($(U("i",i))+$(inds[i])) for i=1:length(inds)]
        if add_index
            return :(C[k,$(subscripts...)])
        else
            return :(C[$(subscripts...)])
        end
    else
        h = symbs[1]
        if length(symbs)>1
            q = symbs[2:end]
        else
            q = []
        end
        exprs = []
        for i = 1:4
            e = Meta.parse( string(h,"_",i,".*(",tensor_prod(q,cat(inds,[i]; dims=1),add_index),")") )
            push!(exprs,e)
        end
        return :(+($(exprs...)))
    end
end

tensor_prod([], Int64[1,2,3,4],)            # C[i1+1,i2+2,i3+3,i4+4]
tensor_prod([], Int64[1,2,3,4], true)      # C[k,i1+1,i2+2,i3+3,i4+4]

tensor_prod(["Phi_1"], Int64[])            # Phi_1_1 * C[i1 + 1] + Phi_1_2 * C[i1 + 2] + Phi_1_3 * C[i1 + 3] + Phi_1_4 * C[i1 + 4]
tensor_prod(["Phi_1", "Phi_2"], Int64[])   # Phi_1_1 * (Phi_2_1 * C[i1 + 1,i2 + 1] + Phi_2_2 * C[i1 + 1,i2 + 2] + Phi_2_3 * C[i1 + 1,i2 + 3] + Phi_2_4 * C[i1 + 1,i2 + 4]) + Phi_1_2 * (Phi_2_1 * C[i1 + 2,i2 + 1] + Phi_2_2 * C[i1 + 2,i2 + 2] + Phi_2_3 * C[i1 + 2,i2 + 3] + Phi_2_4 * C[i1 + 2,i2 + 4]) + Phi_1_3 * (Phi_2_1 * C[i1 + 3,i2 + 1] + Phi_2_2 * C[i1 + 3,i2 + 2] + Phi_2_3 * C[i1 + 3,i2 + 3] + Phi_2_4 * C[i1 + 3,i2 + 4]) + Phi_1_4 * (Phi_2_1 * C[i1 + 4,i2 + 1] + Phi_2_2 * C[i1 + 4,i2 + 2] + Phi_2_3 * C[i1 + 4,i2 + 3] + Phi_2_4 * C[i1 + 4,i2 + 4])
tensor_prod(["Phi_1", "dPhi_2"], Int64[])  # Phi_1_1 * (dPhi_2_1 * C[i1 + 1,i2 + 1] + dPhi_2_2 * C[i1 + 1,i2 + 2] + dPhi_2_3 * C[i1 + 1,i2 + 3] + dPhi_2_4 * C[i1 + 1,i2 + 4]) + Phi_1_2 * (dPhi_2_1 * C[i1 + 2,i2 + 1] + dPhi_2_2 * C[i1 + 2,i2 + 2] + dPhi_2_3 * C[i1 + 2,i2 + 3] + dPhi_2_4 * C[i1 + 2,i2 + 4]) + Phi_1_3 * (dPhi_2_1 * C[i1 + 3,i2 + 1] + dPhi_2_2 * C[i1 + 3,i2 + 2] + dPhi_2_3 * C[i1 + 3,i2 + 3] + dPhi_2_4 * C[i1 + 3,i2 + 4]) + Phi_1_4 * (dPhi_2_1 * C[i1 + 4,i2 + 1] + dPhi_2_2 * C[i1 + 4,i2 + 2] + dPhi_2_3 * C[i1 + 4,i2 + 3] + dPhi_2_4 * C[i1 + 4,i2 + 4])


function create_parameters(d; Tf=Float64)
    lines = []
    for i=1:d
        block = quote
            $(U("M",i)) = orders[$i]
            $(U("start",i)) =  a[$i]
            $(U("dinv",i)) = (orders[$i]-convert($Tf,1.0))/(b[$i]-a[$i])
        end
        for ll in block.args
            push!(lines, ll)
        end
    end
    return lines
end


function create_local_parameters(d; vectorize=true, Tf=Float64)
    lines = []
    if vectorize
        for i=1:d
            bl = quote
                $(U("x",i)) = S[n][$i]
                $(U("u",i)) = ($(U("x",i)) - $(U("start",i)))*$(U("dinv",i))
                $(U("i",i)) = (floor(Int,$(U("u",i)) ))
                $(U("i",i)) = max( min($(U("i",i)),$(U("M",i))-2), 0 )
                $(U("t",i)) = $(U("u",i))-$(U("i",i))
                $(U("tp",i,1)) = $(U("t",i))*$(U("t",i))*$(U("t",i))
                $(U("tp",i,2)) = $(U("t",i))*$(U("t",i))
                $(U("tp",i,3)) = $(U("t",i))
                $(U("tp",i,4)) = 1.0;
            end
            for ll in bl.args
                push!(lines, ll)
            end
        end
    else
        lines = []
        for i=1:d
            bl = quote
                $(U("x",i)) = S[$i]
                $(U("u",i)) = ($(U("x",i)) - $(U("start",i)))*$(U("dinv",i))
                $(U("i",i)) = (floor(Int,$(U("u",i)) ))
                $(U("i",i)) = max( min($(U("i",i)),$(U("M",i))-2), 0 )
                $(U("t",i)) = $(U("u",i))-$(U("i",i))
                $(U("tp",i,1)) = $(U("t",i))*$(U("t",i))*$(U("t",i))
                $(U("tp",i,2)) = $(U("t",i))*$(U("t",i))
                $(U("tp",i,3)) = $(U("t",i))
                $(U("tp",i,4)) = 1.0;
            end
            for ll in bl.args
                push!(lines, ll)
            end
        end
    end
    return lines
end

# function create_local_parameters(d; vectorize=true, Tf=Float64)
#     lines = []
#     if vectorize
#         for i=1:d

#             bl = quote
#                 $(U("x",i)) = S[n][$i]
#                 $(U("u",i)) = ($(U("x",i)) - $(U("start",i)))*$(U("dinv",i))
#                 # $(U("i",i)) = (floor(Int,$(U("u",i)) ))
#                 $(U("i",i)) = unsafe_trunc(Int, $(U("u",i)))
#                 $(U("i",i)) = max( min($(U("i",i)),$(U("M",i))-2), 0 )
#                 $(U("t",i)) = $(U("u",i))-$(U("i",i))
#                 $(U("tp",i,1)) = $(U("t",i))*$(U("t",i))*$(U("t",i))
#                 $(U("tp",i,2)) = $(U("t",i))*$(U("t",i))
#                 $(U("tp",i,3)) = $(U("t",i))
#                 $(U("tp",i,4)) = convert($Tf,1.0);   #TODO
#             end
#             for ll in bl.args
#                 push!(lines, ll)
#             end
#         end
#     else
#         lines = []
#         e = quote
#                 x_ = SVector(S...)
#                 a_ = SVector( a... )
#                 b_ = SVector( b... )
#                 n_ = SVector( orders... )
#                 δ_ = (b_.-a_)./(n_.-1)
#                 # i_ = div.( (x_.-a_), δ_)
#                 i_ = SVector( ()...)
#                 i_ = max.(min.(i_, n_.-2), 0)
#                 λ_ = (x_.-(a_ .+ δ_.*i_))./δ_
#                 ii_ = unsafe_trunc.(Int,i_)

#             end
#         for x in e.args
#             push!(lines, x)
#         end
#         for i=1:d
#             bl = quote


#                 # compat
#                 $(U("i",i)) = ii_[$i]
#                 $(U("t",i)) = λ_[$i]

#                 $(U("tp",i,1)) = $(U("t",i))*$(U("t",i))*$(U("t",i))
#                 $(U("tp",i,2)) = $(U("t",i))*$(U("t",i))
#                 $(U("tp",i,3)) = $(U("t",i))
#                 $(U("tp",i,4)) = convert($Tf,1.0);
#             end
#             for ll in bl.args
#                 push!(lines, ll)
#             end
#         end
#     end
#     return lines
# end

function create_function(d,extrap="natural"; vectorize=true, Tf=Float32)
    if vectorize
        expr = quote
            function $(Symbol(string("eval_UC_spline!")))( a, b, orders, C::Array{T,d}, S::Vector{SVector{d,$Tf}}, V::Vector{T}) where d where T
                $(create_parameters(d; Tf=Tf)...)
                N = size(S,1)
                # @fastmath @inbounds @simd( # doesn't seem to make any difference
                for n=1:N
                    $(create_local_parameters(d;Tf=Tf)...)
                    $(create_Phi(d,extrap,false;Tf=Tf)...)
                    V[n] = $( tensor_prod([string("Phi_",i) for i=1:d], Int64[]; Tf=Tf) )
                end
                # )
            end
        end
    else
        expr = quote
            function $(Symbol(string("eval_UC_spline")))( a, b, orders, C::Array{T,d}, S::SVector{d,U}) where d where T where U
                $(create_parameters(d,Tf=Tf)...)
                N = size(S,1)
                # @fastmath @inbounds @simd( # doesn't seem to make any difference
                $(create_local_parameters(d; vectorize=false,Tf=Tf)...)
                $(create_Phi(d,extrap,false; Tf=Tf)...)
                V = $( tensor_prod([string("Phi_",i) for i=1:d], Int64[]) )
                # )
                return V
            end
        end
    end
    return expr
end
