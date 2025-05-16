
"""
This file is the actual implementation of the information in types.Jl & baseops.Jl
"""

using DelimitedFiles
using FunctionWrappers
import FunctionWrappers: FunctionWrapper
using LinearAlgebra
using Printf
using SparseArrays
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

include("@__DIR__/../type.jl")
include("@__DIR__/../baseops.jl")
include("@__DIR__/../common.jl")
include("@__DIR__/../file_in.jl")
include("@__DIR__/../file_out.jl")
include("@__DIR__/../hc_ham.jl")
include("@__DIR__/../hamiltonian.jl")
include("@__DIR__/../assign.jl")
include("@__DIR__/../ntop.jl")
include("@__DIR__/../transitions.jl")

const csl::Float64 = 29979.2458

BLAS.set_num_threads(Threads.nthreads())

function westereng(molnam::String,ctrl)
# 14 april 25, half the time here is from set up. should look into that
   prm, errs, ℋ, stgs = secordinp(molnam,ctrl)
   ℋ, stgs, errs = opreader(molnam,ctrl,prm,errs,ℋ,stgs)

   σs = σgen_indef(ctrl["NFOLD"])
   σcnt = maximum(size(σs))
   sd = Int(2*ctrl["S"]+1)
   jlist = collect(0.5*iseven(sd):ctrl["Jmax"])
   mcd = Int(2*ctrl["mcalc"]+1)

   #initialize vals
   vtd = ctrl["vtmax"] + 1
   jfd = sd*Int(sum(2.0 .* jlist .+ 1.0))
   vals = zeros(Float64,jfd*vtd,σcnt)
   #initialize vecs
   if ctrl["stages"]==1
      vecs = zeros(Float64,Int(sd*(2*ctrl["Jmax"]+1)*mcd),jfd*vtd,σcnt)
   elseif ctrl["stages"]==2
      vl = sd*(2*ctrl["Jmax"]+1)*(ctrl["mmax"]+1)
      vecs = zeros(Float64,Int(vl),jfd*vtd,σcnt)
   end
   #initialize tvecs
   tvecs = zeros(2*ctrl["mcalc"]+1,ctrl["mmax"]+1,σcnt)
#   #do the energy level calculation
   if ctrl["stages"]==1
      @time vals,vecs = tsrcalc_1stg!(vals,vecs,jlist,σs,ctrl,prm,stgs,ℋ)
   elseif ctrl["stages"]==2
      tvals = zeros(ctrl["mmax"]+1,σcnt)
@time vals,vecs,tvals,tvecs = tsrcalc_2stg!(vals,vecs,tvals,tvecs,jlist,σs,ctrl,prm,stgs,ℋ)
   else
      @warn "Invalid stages number"
   end
   qns = bigqngen(size(vals,1),jlist,ctrl["S"],ctrl["vtmax"],σs)
   #outputs!
   if occursin("E",ctrl["RUNmode"])
      engwriter(molnam,ctrl,vals,qns)
   end
   if iszero(tvecs)
      tvecs = [0.0]
   end
   return vals,vecs,qns, tvecs
end

function westersim(molnam,ctrl,fvls,fvcs,fqns,tvecs)
   μf = intreader(molnam,ctrl) # function doesn't exist yet
   σcnt = σcount(ctrl["NFOLD"])
   jmax = ctrl["Jmax"]
   #perm = permdeterm(scls,stg)
   kbT = ctrl["TK"]*20836.61912 #MHz/K
   Qrt = qrtcalc(fvls,ctrl["TK"])
   finfrq = zeros(0,4)
   finqns = zeros(Int,0,12)
   for sc ∈ 1:σcnt
      σ = sc - 1
      vals = fvls[:,sc]
      vecs = fvcs[:,:,sc]
      quns = fqns[:,:,sc]
      if ctrl["stages"]==1
         #add uncertainty calculator
         fr, qn = tracalc(ctrl,vals,vecs,tvecs,quns,μf,σ,Qrt)
      elseif ctrl["stages"]==2
         #add uncertainty calculator
         @warn "function under construction"
      else
         @warn "this number of stages isn't implemented"
      end
      finfrq = vcat(finfrq,fr)
      finqns = vcat(finqns,qn)
   end#sigma loop
   #TraWriterSPCAT(molnam, ctrl["S"], finfrq, finqns)
   #TraWriter(molnam, ctrl["S"], finfrq, finqns)
   return finfrq, finqns
end#westersim


function westermain(molnam::String)
   molnam = String(split(molnam,'.')[1])
   #read input file
   ctrl = ctrlinp(molnam)
   if occursin("T",ctrl["RUNmode"])
      println("This RUNmode is probably not what you want")
      westersim(molnam,nothing, ctrl)
      westerfit(molnam, ctrl)
   else
      if occursin("F",ctrl["RUNmode"])
         prm, pcov = westerfit(molnam, ctrl)
      else
         prm = nothing
         puncs = nothing
         pcov = nothing
      end
      if occursin("E", ctrl["RUNmode"])||occursin("S", ctrl["RUNmode"])
         vas,ves,qns,μs,prm,scls,stg,cdo,tvcs = westereng(molnam, prm, ctrl)
         if occursin("S", ctrl["RUNmode"])
            westersim(molnam,prm,ctrl,vas,ves,qns,μs,prm,scls,stg,cdo,pcov,tvcs)
         end
      end
   end
end

#Base.@ccallable function main()::Cint
function main(molnam)
   #molanm = ARGS[1]
   ctrl = blockfind_all(molnam)
   ctrl = ctrlinp(molnam,ctrl)
   vals,vecs,qns,tvecs = westereng(molnam,ctrl)
   westersim(molnam,ctrl,vals,vecs,qns,tvecs)
   return 0
end
