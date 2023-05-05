using DelimitedFiles

function reorder()
line = readdlm("temp2.lne",',')
perma = [4, 1, 2, 3, 8, 5, 6, 7, 9, 10]
perme = [4, 1, 2, 3, 8, 5, 6, 7, 11, 12]
linea = line[:,perma]
linee = line[:,perme]
linetotal = vcat(linea,linee)
linetotal[:,1] = linetotal[:,1] ./ 2
linetotal[:,5] = linetotal[:,5] ./ 2
writedlm("ocltol.lne", linetotal, ',')
end

reorder()