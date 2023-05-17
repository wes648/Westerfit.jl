using DelimitedFiles

function reorder()
line = readdlm("ipfp.temp",',')
line = hcat(line,zeros(Int,size(line,1),2))
perma = [4, 1, 2, 3, 11, 8, 5, 6, 7, 12, 9, 10]
linetotal = line[:,perma]
linetotal[:,1] = linetotal[:,1] .- 0.5
linetotal[:,6] = linetotal[:,6] .- 0.5
writedlm("ipfp.lne", linetotal, ',')
end

reorder()
