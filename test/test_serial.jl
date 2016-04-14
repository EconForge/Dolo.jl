N = 3
mats_a = [rand(2,3) for n=1:N]
mats_b = [rand(3,4) for n=1:N]
mats_c = [mats_a[i]*mats_b[i] for i=1:N]

A = cat(1,[reshape(a,1,2,3) for a in mats_a]...)
B = cat(1,[reshape(a,1,3,4) for a in mats_b]...)
C = serial_multiply(A,B)

diff = C - cat(1,[reshape(a,1,2,4) for a in mats_c]...)

@test maximum(abs(diff))==0
