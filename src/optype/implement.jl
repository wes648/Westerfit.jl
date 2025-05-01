
"""
This file is the actual implementation of the information in types.Jl & baseops.Jl
"""


using DelimitedFiles
using FunctionWrappers
import FunctionWrappers: FunctionWrapper
using LinearAlgebra
using Printf
using SparseArrays
using WIGXJPFjl
#using JET
#using BenchmarkTools
using ProfileView

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


function westereng(molnam::String)
# 14 april 25, half the time here is from set up. should look into that
   ctrl = blockfind_all(molnam)
   ctrl = ctrlinp(molnam,ctrl)
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

   μf = [μOb(1.0,μzf,Et_int,1,0,0)]
   intcalc(ctrl,vecs[:,:,1],μf,0)
   return vals,vecs,qns, tvecs
end

function westersim(vals,vecs,qns,tvecs)
   μs = intreader(molanm,ctrl) # function doesn't exist yet
end

