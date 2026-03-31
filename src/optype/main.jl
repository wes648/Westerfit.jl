
"""
This file is the actual implementation of the information in types.jl & baseops.jl
"""

using DelimitedFiles
using FunctionWrappers
import FunctionWrappers: FunctionWrapper
using LinearAlgebra
import Printf: @sprintf
#import speed of light, plank, atomic mass, electron mass, fine structure constant
import PhysicalConstants.CODATA2022: c_0
using SparseArrays
using TOML
@static if Sys.iswindows()
   using WignerSymbols
   wig3j(a,b,c,d,e,f) = wigner3j(Float64,a,b,c,d,e,f)
   wig6j(a,b,c,d,e,f) = wigner6j(Float64,a,b,c,d,e,f)
else
   using WIGXJPFjl
end
#using JET
#using BenchmarkTools
#using ProfileView
include("@__DIR__/../file_in.jl")
include("@__DIR__/../psi.jl")
include("@__DIR__/../type.jl")
include("@__DIR__/../baseops.jl")
include("@__DIR__/../common.jl")
#include("@__DIR__/../derivatives.jl")
#include("@__DIR__/../file_out.jl")
#include("@__DIR__/../hc_ham.jl")
include("@__DIR__/../hamil.jl")
include("@__DIR__/../assign.jl")
include("@__DIR__/../ntop.jl")
#include("@__DIR__/../transitions.jl")
#include("@__DIR__/../opt-com.jl")
#include("@__DIR__/../optimizer.jl")

#const csl::Float64 = 29979.2458
const csl::Float64 = (c_0 * 1e-4).val #MHz/cm⁻¹
#csl = 29979.2458
const hccount::Int = 11
BLAS.set_num_threads(Int(0.5*Sys.CPU_THREADS))
@show Threads.nthreads()

function westereng(ctrl,prm,ℋ)::Eigs 
   wvs = Eigs(ctrl)
   H_calc(ctrl,wvs,prm,ℋ)
   return wvs
end

function westermain()
   molnam = "test"
   ctrl, ℋ, prms, scls = inp_reader(molnam)
   wvs = westereng(ctrl,prms, ℋ)
   return wvs
end

function westersim(molnam,ctrl,fvls,fvcs,fqns,tvecs)
   μf = intreader(molnam,ctrl) # function doesn't exist yet
   σcnt = σcount(ctrl.NFOLD)
   jmax = ctrl.Jmax
   #perm = permdeterm(scls,stg)
   kbT = ctrl.TK*20836.61912 #MHz/K
   Qrt = qrtcalc(fvls,ctrl.TK)
   finfrq = zeros(0,4)
   finqns = zeros(Int,0,12)
   for sc ∈ 1:σcnt
      σ = sc - 1
      vals = fvls[:,sc]
      vecs = fvcs[:,:,sc]
      quns = fqns[:,:,sc]
      if ctrl.stages==1
         #add uncertainty calculator
         fr, qn = tracalc(ctrl,vals,vecs,tvecs,quns,μf,σ,Qrt)
      elseif ctrl.stages==2
         #add uncertainty calculator
         @warn "function under construction"
      else
         @warn "this number of stages isn't implemented"
      end
      finfrq = vcat(finfrq,fr)
      finqns = vcat(finqns,qn)
   end#sigma loop
   #TraWriterSPCAT(molnam, ctrl.S, finfrq, finqns)
   #TraWriter(molnam, ctrl.S, finfrq, finqns)
   return finfrq, finqns
end#westersim


function westermain_old(molnam::String)
   molnam = String(split(molnam,'.')[1])
   #read input file
   ctrl = ctrlinp(molnam)
   if occursin("T",ctrl.RUNmode)
      println("This RUNmode is probably not what you want")
      westersim(molnam,nothing, ctrl)
      westerfit(molnam, ctrl)
   else
      if occursin("F",ctrl.RUNmode)
         prm, pcov = westerfit(molnam, ctrl)
      else
         prm = nothing
         puncs = nothing
         pcov = nothing
      end
      if occursin("E", ctrl.RUNmode)||occursin("S", ctrl.RUNmode)
         vas,ves,qns,μs,prm,scls,stg,cdo,tvcs = westereng(molnam, prm, ctrl)
         if occursin("S", ctrl.RUNmode)
            westersim(molnam,prm,ctrl,vas,ves,qns,μs,prm,scls,stg,cdo,pcov,tvcs)
         end
      end
   end
end


function westerfit(molnam::String,ctrl::Controls)
"""
   The fitter!
"""
   println("westerfit!")
   prm, errs, stgs = secordinp(molnam,ctrl)
   ℋ, stgs, errs, unts = opreader(molnam,ctrl,prm,errs,stgs)
   #if occursin("F",ctrl.RUNmode) #Normal Fit behavior, overpowers check
      lines = readdlm("$molnam.lne", ',', Float64,comments=true,comment_char='#',
                     skipblanks=true)
      linelength = (size(lines,1))
#   else # Self-consistency check
#      lines = readdlm("$molnam.cat", ',', Float64,comments=true,comment_char='#')
#      lines = pred2lne(lines,ctrl.S)
#   end

   #determine the states
   linds, ofreqs, luncs = lineprep(lines,ctrl.NFOLD[1],ctrl.S,ctrl.vtmax)
   #@show linds
   jlist = jlister(linds)
   #opt
   #outputinit(molnam,prm,errs,linelength,ctrl)
   outputinit2(molnam,prm,errs,linelength,ctrl,ℋ)#16.0,8.0)
   tsrp, pcov, omcs, cfrqs, vals = lbmq(ctrl,jlist,ofreqs,luncs,linds,
                                             prm,errs,ℋ,stgs,molnam)
   reswritter(molnam,lines,omcs,cfrqs)
   return tsrp, pcov
end


#Base.@ccallable function main()::Cint
function main(molnam)
   #molanm = ARGS[1]
   @time ctrl = blockfind_all(molnam) #0.3 s
   @time ctrl = ctrlinp(molnam,ctrl) # 1s
   #vals,vecs,qns,tvecs = westereng(molnam,ctrl)
   #westersim(molnam,ctrl,vals,vecs,qns,tvecs)
   westerfit(molnam,ctrl)
   return 0
end
