
using DelimitedFiles
include("@__DIR__/../../../src/main.jl")


const csl = 29979.2458
qns = ["N"; "Ka"; "Kc"; "m"; "σ"; "E"]

bash(str::String) = run(`bash -c $str`)
macro sh_str(str::String)
   bash(str)
end

function setvars(nams,vals)
   run(`cp belgi.ref input.txt`)
   run(`cp ref.inp belgi.inp`)
   wvls = copy(vals)
   wvls[1:6] *= csl
   for i in 1:length(vals)
      nam = nams[i]
      val = vals[i]
      wvl = wvls[i]
      sedstr = "s/$nam/$val/g"
      belg = `sed -i $sedstr input.txt`
      #bstr = "sed -i 's/$nam/$val/g' input.txt"
      run(belg,wait=true)
      sedstr = "s/$nam/$wvl/g"
      west = `sed -i $sedstr belgi.inp`
      run(west,wait=true)
   end
end

#run BELGI-CS 
function runbelgi()
   run(pipeline(`wine $(homedir())/rot/bin/belgi-cs`, stdout="output"),wait=true)
   #run(`cp output belgi.output`)
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
   run(`sed -i 's/*//g' output`)
   run(`sed -i 's/PURITY=/,  1, 1/g' output`, wait=true)
end

function vt2m(vt,σ)
   if vt==zero(vt)
      return σ
   elseif iseven(vt)
      return 3*(vt-1) + σ
   elseif isodd(vt)
      return -3*vt + σ
   end
end

#process the belgi output
function procbelgi()
   belg = readdlm("output",',')
   #re-order to resemble westerfit
   #       N, K, k, m, σ, E
   perm = [4, 2, 5, 1, 6, 3]
   belg = belg[:,perm]
   #println(belg)
   #convert vt to m
   @. belg[:,3] = vt2m(belg[:,3], belg[:,5])
   #sort by energy
   belg = belg[sortperm(belg[:,end]), :]
   return belg
end

function rms(a::Array,b::Array)::Float64
   c = a .- b
   return BLAS.nrm2(c) / √(length(c))
end

function westvbelg(belg,west)
   west = west[sortperm(west[:,end]),:]
   for i in 1:6
      err = rms(belg[:,i], west[:,i])
      println("RMS of $(qns[i]) = $err")
   end
end

function runtest()
   nams = ["Av"; "Bv"; "Cv"; "Dabv"; "Fv"; "V3v"; "rhov"]
   vals = [3000.0/csl; 1500.0/csl; 1000.0/csl; 0.0/csl; 5.1; 200.0; 0.05]
   setvars(nams,vals)
   runbelgi()
   westerfit("belgi")
   belg = procbelgi()
   west = readdlm("belgi.eng",',')
   westvbelg(belg,west)
end
#

runtest()
#ΔE = sort(belg[:,3]) .- sort(west[:,6])
