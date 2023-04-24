#prevare var file
#cp blank.var test.var
##copy template file
##set spin degeneracy
#sed -i 's/SD/$SD/g' test.var
##set K max
##replace varameter values
#prevare int file
#cp blank.int test.inp
##set temperature
##set Fmax
##set νmax
#run simulation
#grab vartion function value
#replace vartition function value
#rerun simulation
#convert cat file to csv
##should be similar to my pred2csv script
#convert csv to lne
#build westerfit varameter input
#fit spcat output in westerfit
#print results summary
#rm test.*

using DelimitedFiles
using LinearAlgebra

qns = ["J", "N", "Ka", "Kc", "E"]

function setvars(nams,vals,s) #sets all the variables to the values in nams and vals
	run(`cp spfit.ref test.var`)
	run(`cp ref.inp spfit.inp`)
	run(`sed -i "s/sv/$s/g" spfit.inp`)
	sd = Int(2s+1)
	run(`sed -i "s/sd/$sd/g" test.var`)
	for i in 1:length(vals)
		nam = nams[i]
		val = vals[i]
		sedstr = "s/$nam/$val/g"
		west = `sed -i $sedstr spfit.inp`
		run(west)
		if nam=="chizz"
			val *= 1.5
		elseif nam=="chixmy"
			val *= 0.25
		end
		sedstr = "s/$nam/$val/g"
		sp = `sed -i $sedstr test.var`
		run(sp)
	end
end

#run spfit
function runspfit()
	run(`spcat test.var`) #runs spcat
	#run(`gawk -i inplace '$1=$1' FIELDWIDTHS='13 8 8 2 10 3 7 6 2 4 6 2' OFS=, test.cat`)
	run(`sed -i "s/://g" test.egy`) #gets rid of any colons in test.egy
	run(`gawk -i inplace '$1=$1' FIELDWIDTHS='6 5 18 18 11 5 3 3 3' OFS=, test.egy`) #comma delimited
end

function runwesterfit() #runs westerfit
	run(`julia -t4 ../../src/run.jl spfit`) #runs westerfit
end

function procspfit()
	spfit = readdlm("test.egy",',')
	perm = [6, 7, 8, 9, 3]
	spfit = spfit[:,perm]
	spfit = spfit[sortperm(spfit[:,end]), :]
	spfit[:,1] = (spfit[:,1].-1)./2
	return spfit
end

function procwesterfit()
	west = readdlm("spfit.eng", ',')
	perm = [1, 2, 3, 4, 7]
	west = west[:,perm]
	west = west[sortperm(west[:,end]), :]
	west[:,1] = west[:,1].*0.5
	west[:,3] = abs.(west[:,3])
	return west
end

function rms(a::Array,b::Array)::Float64
   c = a .- b
   return BLAS.nrm2(c) / √(length(c))
end

function westvspfit(spfit,west)
	# J N Ka Kc E
	for i in 1:5
		err = rms(spfit[:,i], west[:,i])
		println("RMS of $(qns[i]) = $err")
	end
end

function runtest()
	nams = ["Av"; "Bv"; "Cv"; "DK"; "DNK"; "DN"; "dN"; "dK"; "chizz"; "chixmy"; "chixz"]
	vals = [3000.0; 1500.0; 1000.0; .3; .15; .2; .1; .08; 20.; 10.; 5.]
	s = 1.0
	setvars(nams, vals, s)
	runspfit()
	runwesterfit()
	west = procwesterfit()
	spfit = procspfit()
	westvspfit(spfit,west)
end

runtest()