using AxisArrays
using Unitful
import Unitful: s, ms, µs

fs = 40000 # Generate a 40kHz noisy signal, with spike-like stuff added for testing
y = randn(60*fs+1)*3
for spk = (sin.(0.8:0.2:8.6) .* [0:0.01:.1; .15:.1:.95; 1:-.05:.05]   .* 50,
         sin.(0.8:0.4:8.6) .* [0:0.02:.1; .15:.1:1; 1:-.2:.1] .* 50)
  i = rand(round(Int,.001fs):1fs)
  while i+length(spk)-1 < length(y)
      y[i:i+length(spk)-1] += spk
      i += rand(round(Int,.001fs):1fs)
  end
end

y
[y 2y]

A = AxisArray([y 2y], Axis{:time}(0s:1s/fs:60s), Axis{:chan}([:c1, :c2]))

a =A[Axis{:time}(4)]

collect(a)
