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
using Plots
using LinearAlgebra
include("@__DIR__/../../../src/main.jl")

qns = ["J", "N", "Ka", "Kc", "E"]
csl = 29979.2458

function setvars(nams,vals,s) #sets all the variables to the values in nams and vals
   run(`cp ref.var test.var`)
   run(`cp ref.inp test.inp`)
   run(`sed -i "s/sv/$s/g" test.inp`)
   sd = Int(2s+1)
   run(`sed -i "s/sd/$sd/g" test.var`)
   for i in 1:length(vals)
      nam = nams[i]
      val = vals[i]
      sedstr = "s/$nam/$val/g"
      west = `sed -i $sedstr test.inp`
      run(west)
      if nam=="chizz"
         val *= 1.5
      elseif nam=="chixmy"
         val *= 0.25
      elseif nam=="chixz"
         val *= -1.0
      #elseif nam=="epxz"
      #   val *= -1.0
      #else
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
   #run(`julia -t4 ../../src/run.jl spfit`) #runs westerfit
   westerfit("test")
end

function procspfit()
   spfit = readdlm("test.egy",',')
   perm = [6, 7, 8, 9, 3]
   spfit = spfit[:,perm]
   #spfit = spfit[sortperm(spfit[:,end], by=abs), :]
#   spfit = spfit[sortperm(spfit[:,1]), :]
   spfit = spfit[sortperm(spfit[:,5]), :]
   #spfit = spfit[sortperm(spfit[:,1]), :]
   spfit[:,1] = (spfit[:,1].-1)./2
   spfit[:,end] ./= csl
   return spfit
end

function procwesterfit()
   west = readdlm("test.eng", ',')
   perm = [1, 2, 3, 4, 7]
   west = west[:,perm]
   #west = west[sortperm(west[:,end], by=abs), :]
   west = west[sortperm(west[:,5]), :]
   #west = west[sortperm(west[:,1]), :]
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
   err = rms(spfit[:,end], west[:,end])*csl
   println("RMS of $(qns[end]) = $err")
end

function runtest()
   nams = ["Av"; "Bv"; "Cv"; "Dabv"; "chizz"; "chixmy"; "chixz"; "epzz"; "epxx"; "epyy"; "epxz"]
   #vals = [3000.0; 1500.0; 1000.0; 50.0; -30.; 20.; 50.; 000.; 00.; 0.; 00.] #just qua
   vals = [3000.0; 1500.0; 1000.0; 50.0; -30.; -20.; 50.; 300.; 80.; 13.; 100.] #all
   #vals = [3000.0; 1500.0; 1000.0; 50.0; -00.; 00.; 00.; 300.; 80.; 13.; 100.] #just sr
   #vals = [3000.0; 1500.0; 1000.0; 0.0; -30.; 20.; 0.; 300.; 80.; 13.; 00.] #on-diag
   s = 1.0
   setvars(nams, vals, s)
   runspfit()
   runwesterfit()
   west = procwesterfit()
   spfit = procspfit()
   westvspfit(spfit,west)
#   println([west[:,1] spfit[:,1]])
   plot(west[:,5],(spfit[:,5].-west[:,5]).*csl,seriestype=:scatter,label=false,ylab="ΔE (MHz)",
      xlab="westerfit Energy (cm⁻¹)",dpi=500)
   hline!([-1e-6,1e-6],label=false)
   savefig("wesvspft.png")
end

runtest()

#=
function stabdiv(a,b)
   if b==zero(b)
      return 0.0
   else
      return a ./ b
   end
end

function atest(χa)
   nams = ["Av"; "Bv"; "Cv"; "Dabv"; "DK"; "DNK"; "DN"; "dN"; "dK"; "chizz"; "chixmy"; "chixz"; "epzz"; "epxx"; "epyy"; "epxz"]
   #vals = [3000.0; 1500.0; 1000.0; .3; .15; .2; .1; .08; 20.; 10.; 5.]
   #vals = [3000.0; 1500.0; 1000.0; 00.0; .0; .0; .0; .0; .0; χa; 400.; 15.]
   vals = [000.0; 00.0; 000.0; 0.0; .0; .0; .0; .0; .0; 0.; 0.; 0.; χa; 00.; 0.; 00.]
   s = 1.0
   setvars(nams, vals, s)
   runspfit()
   runwesterfit()
   west = procwesterfit()
   spfit = procspfit()
   westvspfit(spfit,west)
   omc = spfit[:,end] .- west[:,end]
   println("when ϵaa = $χa, sum(omc) = $(sum(omc))")
   omc = stabdiv(omc, χa)
   return omc
end

function runeatest()
   eas = 15.0 .* (10.0 .^collect(-2:2))
   eas = vcat([0.0], eas)
   len = length(eas)
   omcs = zeros(Float64,75,len) 
   for i in 1:len
      omcs[:,i] .= atest(eas[i])
   end
   x = collect(1:75)
   plot(x, omcs[:,1], seriestype=:scatter)
   for i in 2:len
      plot!(x, omcs[:,i], seriestype=:scatter)
   end
   savefig("parm.png")
end
=#