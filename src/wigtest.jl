
include("./WIGXJPF.jl")
using .WIGXJPF

out = zeros(Float64,10)

Threads.@threads for i in 1:10
   out[i] = WIGXJPF.wig3j(10,9,i,i,0,-i)
end

println(out)
