
using DelimitedFiles
using FunctionWrappers
import FunctionWrappers: FunctionWrapper
using LinearAlgebra
import Printf: @sprintf
#import speed of light, plank, atomic mass, electron mass, fine structure constant
import PhysicalConstants.CODATA2022: c_0, k_B, h
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
include("@__DIR__/../common.jl")
#include("@__DIR__/../derivatives.jl")
include("@__DIR__/../file_out.jl")
include("@__DIR__/../hc_ham.jl")
include("@__DIR__/../hamil.jl")
include("@__DIR__/../assign.jl")
include("@__DIR__/../ntop.jl")
include("@__DIR__/../transitions.jl")
#include("@__DIR__/../opt-com.jl")
#include("@__DIR__/../optimizer.jl")

include("@__DIR__/../baseops.jl")
include("@__DIR__/../dipoles.jl")

#const csl::Float64 = 29979.2458
const csl::Float64 = (c_0 * 1e-4).val # MHz / cm⁻¹
const kb::Float64 = (k_B/h * 1e-6).val # MHz / K
#csl = 29979.2458
const hccount::Int = 11
BLAS.set_num_threads(Int(0.5*Sys.CPU_THREADS))
@show Threads.nthreads()

function westereng(molnam, ctrl,prm,ℋ)::Eigs 
   wvs = Eigs(ctrl)
   jsσs = jσlister_full(ctrl.S,ctrl.Jmax, σcount(ctrl.NFOLD))
   H_calc(ctrl,wvs,prm,ℋ,jsσs)
   if occursin("E",ctrl.RUNmode)
      engwriter(molnam, ctrl.Jmax, ctrl.S, ctrl.vtmax, wvs.rst.vals)
      println("yay! energy levels writen to $molnam.eng")
   end
   return wvs
end

function westersim(molnam, ctrl, μs, wvs)
   σs = σcount(ctrl.NFOLD)
   #jlst = jbjk
   frqs, inds = tracalc(ctrl,μs,wvs)
#   @show size(wvs.rst.vals)
#   @show inds
   writefreqs(molnam,ctrl,frqs,inds)
end

function westermain()
   molnam = "test"
   @time ctrl, ℋ, prms, scls, μs = inp_reader(molnam)
   @time wvs = westereng(molnam, ctrl,prms, ℋ)
   @time westersim(molnam, ctrl, μs, wvs)
   return wvs
end



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
