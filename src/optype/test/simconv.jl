
using DelimitedFiles
sim = readdlm("./old_rt/2ba.sim",',')
l = size(sim[:,1])
writedlm("test_rt.lne", [sim[:,1] sim[:,1:3] zeros(l) sim[:,4:5] sim[:,5:7] zeros(l) sim[:,8:9] fill(0.000001,l) ], ',')
