
using Base.Threads
include("./WIGXJPF.jl")
using .WIGXJPF

nthr = nthreads()
out = zeros(Float64,nthr)

println("Serial")
for i in 1:nthr
   out[i] = wig3j(nthr,2,nthr,i,0,-i)
end
println(out)


println("Parallel")
out = zeros(Float64,nthr)
@threads for i in 1:nthr
   out[i] = wig3j(nthr,2,nthr,i,0,-i)
end
println(out)
