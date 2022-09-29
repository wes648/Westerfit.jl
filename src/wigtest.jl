
using Base.Threads
include("./WIGXJPF.jl")
using .WIGXJPF
#using WignerSymbols

nthr = nthreads()

out = zeros(Float64,nthr)

@threads for i in 1:nthr
   out[i] = wig3j(Float64,nthr,2,nthr,i,0,-i)
end

println(out)
