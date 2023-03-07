
using DelimitedFiles

const csl = 29979.2458
nams = ["Av"; "Bv"; "Cv"; "Dabv"; "Fv"; "V3v"; "ρv"]
vals = [3000.0; 1500.0; 1000.0; 50.0; 5.0*csl; 200*csl; 0.05]

function setvars(strs,vals)
	run(`cp belgi.ref input.txt`)
	for i in 1:length(vals)
		belg = `sed -i 's/$(strs[i])/$(vals[i])/g input.txt`
		run(belg)
		west = `sed -i 's/$(strs[i])/$(vals[i])/g belg.inp`)
		run(wes)
	end
end

#run BELGI-CS 
function runbelgi()
	run(pipeline(`wine /home/wes/rot/bin/belgi-cs`, stdout="output"),wait=true)
	#EXTRACT ENERGY LEVELS
	run(`sed -i -n "/BEGIN DO/,/END/p" output`, wait=true)
	#FILTER OUT OTHER LINES
	run(`sed -i '/CM-1/!d' output`)
	run(`sed -i 's/VT=//g' output`)
	run(`sed -i 's/K=/,/g' output`)
	run(`sed -i 's/A=/,/g' output`)
	run(`sed -i 's/E=/,/g' output`)
	run(`sed -i 's/CM-1  N=/,/g' output`)
	run(`sed -i 's/PAR= +/,  1, 0/g' output`)
	run(`sed -i 's/PAR= -/, -1, 0/g' output`)
	run(`sed -i 's/*//,  1, 1' output`)
	run(`sed -i 's/PURITY=/, 1/g' output`)
end


#process the belgi output
function procbelgi()
	belg = readdlm("output",',')
	#re-order to resemble westerfit
	perm = [4, 2, 5, 1, 6, 3]
	belg = belg[perm,:]
	#convert vt to m

	#sort by energy
	perm = sortperm(belg[:,end])
	belg = belg[:,perm]
	return belg
end
#west = readdlm("belgi.eng",',')

#ΔE = sort(belg[:,3]) .- sort(west[:,6])
