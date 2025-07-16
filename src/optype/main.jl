
"""
This file is the actual implementation of the information in types.Jl & baseops.Jl
"""

using DelimitedFiles
using FunctionWrappers
import FunctionWrappers: FunctionWrapper
using LinearAlgebra
import Printf: @sprintf
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
include("@__DIR__/../derivatives.jl")
include("@__DIR__/../file_in.jl")
include("@__DIR__/../file_out.jl")
include("@__DIR__/../hc_ham.jl")
include("@__DIR__/../hamiltonian.jl")
include("@__DIR__/../assign.jl")
include("@__DIR__/../ntop.jl")
include("@__DIR__/../transitions.jl")
include("@__DIR__/../opt-com.jl")
include("@__DIR__/../optimizer.jl")

const csl::Float64 = 29979.2458

BLAS.set_num_threads(Threads.nthreads())
@show Threads.nthreads()

function westereng(molnam::String,ctrl::Controls)
# 14 april 25, half the time here is from set up. should look into that
   prm, errs, ℋ, stgs = secordinp(molnam,ctrl)
   ℋ, stgs, errs = opreader(molnam,ctrl,prm,errs,ℋ,stgs)

   σs = σgen_indef(ctrl.NFOLD)
#   @show σs
   σcnt = maximum(size(σs))
#   @show σcnt
   sd = Int(2*ctrl.S+1)
#   jlist = collect(0.5*iseven(sd):ctrl.Jmax)
   jlist = collect(iseven(sd):2:Int(2*ctrl.Jmax))
   jlist = [kron(jlist,ones(Int,σcnt)) kron(ones(Int,length(jlist)),collect(1:σcnt))]
   mcd = Int(2*ctrl.mcalc+1)

   #initialize vals
   vtd = ctrl.vtmax + 1
#   @show sd
#   @show jlist
   jfd = sd*Int(sum(jlist[:,1] .+ 1.0))
#   @show jfd
   vals = zeros(Float64,jfd*vtd,σcnt)
   #initialize vecs
   if ctrl.stages==1
      tvecs = zeros(1)
      vecs = zeros(Float64,Int(sd*(2*ctrl.Jmax+1)*mcd),jfd*vtd,σcnt)
   elseif ctrl.stages==2
      vl = sd*(2*ctrl.Jmax+1)*(ctrl.mmax+1)
      vecs = zeros(Float64,Int(vl),jfd*vtd,σcnt)
      #initialize tvecs
      tvecs = zeros(2*ctrl.mcalc+1,ctrl.mmax+1,σcnt)
   end
#   #do the energy level calculation
   if ctrl.stages==1
      @time vals,vecs = tsrcalc_1stg!(vals,vecs,jlist,σs,ctrl,prm,stgs,ℋ)
   elseif ctrl.stages==2
      tvals = zeros(ctrl.mmax+1,σcnt)
@time vals,vecs,tvals,tvecs = tsrcalc_2stg!(vals,vecs,tvals,tvecs,jlist,σs,ctrl,prm,stgs,ℋ)
   else
      @warn "Invalid stages number"
   end
   qns = bigqngen(size(vals,1),jlist,ctrl.S,ctrl.vtmax,σs)
   #outputs!
#   @show size(vals)
   if occursin("E",ctrl.RUNmode)
      engwriter(molnam,ctrl,vals,qns)
   end
   if iszero(tvecs)
      tvecs = [0.0]
   end
   return vals,vecs,qns, tvecs
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


function westermain(molnam::String)
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
   prm, errs, ℋ, stgs = secordinp(molnam,ctrl)
   ℋ, stgs, errs, unts = opreader(molnam,ctrl,prm,errs,ℋ,stgs)
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
   tsrp, pcov, omcs, cfrqs, vals = lbmq(ctrl,jlist,ofreqs,luncs,linds,
                                             prm,errs,ℋ,stgs,molnam)
   reswritter(molnam,lines,omcs,cfrqs)
   return tsrp, pcov
end


#Base.@ccallable function main()::Cint
function main(molnam)
   #molanm = ARGS[1]
   @time ctrl = blockfind_all(molnam)
   @time ctrl = ctrlinp(molnam,ctrl)
   #vals,vecs,qns,tvecs = westereng(molnam,ctrl)
   #westersim(molnam,ctrl,vals,vecs,qns,tvecs)
   westerfit(molnam,ctrl)
   return 0
end
