
using MKL
using DelimitedFiles
using LinearAlgebra
using LinearAlgebra.BLAS
using Printf
using SparseArrays
using StaticArrays
using Base.Threads
using Dates
#include("@__DIR__/../WIGXJPF.jl")

@static if Sys.iswindows()
   using WignerSymbols
   wig3j(a,b,c,d,e,f) = wigner3j(Float64,a,b,c,d,e,f)
   wig6j(a,b,c,d,e,f) = wigner6j(Float64,a,b,c,d,e,f)
else
   using WIGXJPFjl
end

include("@__DIR__/../assign.jl")
include("@__DIR__/../common.jl")
include("@__DIR__/../files_in.jl")
include("@__DIR__/../files_out.jl")
#include("@__DIR__/../jacobi.jl")
include("@__DIR__/../multstage.jl")
include("@__DIR__/../new_ham.jl")
include("@__DIR__/../lsq/opt-com.jl")
include("@__DIR__/../lsq/optimizer.jl")
include("@__DIR__/../lsq/opt-approx.jl")
include("@__DIR__/../lsq/opt-dubious.jl")
#include("@__DIR__/../opt/opt-twostg.jl")
include("@__DIR__/../remnant.jl")
include("@__DIR__/../transitions.jl")
include("@__DIR__/../transi-2stg.jl")
#include("@__DIR__/../userdefop.jl")

BLAS.set_num_threads(Threads.nthreads())

const csl = 29979.2458

function westereng(molnam::String, prm, ctrl)
   println("westersim!")
   molnam = String(split(molnam,'.')[1])
   #molnam = replace(molnam, r".inp"=>"")
   #read input file
   sof, ser = secordinp(molnam,String(ctrl["Irrep"]))
   μs, cdf, cdn, cde, cdo, stg = opinp(molnam,ctrl["NFOLD"])
   if prm==nothing
      prm = vcat(sof,cdf)
   else
      #this occurs when westerfit is called before thus feeding the new
      # parameters in to the function in their natural form
      #the parameters aren't being double transformed. it's fine i swear
      prm[1:18] = sod2prep_full(prm[1:18])
   end
   #calculate energy levels
   σcnt = σcount(ctrl["NFOLD"])
   sd = Int(2.0*ctrl["S"]+1.0)
   jmin = 0.5*iseven(sd)
   jmax = ctrl["Jmax"]
   jlist = collect(Float64,jmin:jmax)
   vtd = Int(ctrl["vtmax"]+1)
   if ctrl["stages"] == 1
      mcd = Int(2*ctrl["mcalc"]+1)
      tvcs = Array{Float64,3}(undef,0,0,0)
   elseif ctrl["stages"] == 2
      mcd = Int(ctrl["mmax"]+1)
      tvcs = zeros(2*ctrl["mcalc"]+1,mcd,σcnt)
   else
      println("DANGER")
   end

   jfd = sd*Int(sum(2.0 .* jlist .+ 1.0))
   fvls = zeros(Float64,jfd*vtd,σcnt)
   #mcd = Int(2*ctrl["mcalc"]+1)
   fqns = zeros(Int,jfd*vtd,6,σcnt)
   fvcs = zeros(Float64,Int(sd*(2*jmax+2)*mcd),jfd*vtd,σcnt)
   @time for sc in 1:σcnt
      σ = sc - 1
      jmsd = Int(mcd*sd*(2*jmax+1))
      jsvd = Int(jfd*vtd)
#       fvls[1:jvsd,σ], fvcs[1:jmsd,1:jsvd,σ], fqns[:,1:jvsd,σ] = tsrcalc(sof,cdf,cdo,
   if ctrl["stages"] == 1
      tempa, tempe, tempq = tsrcalc(ctrl,prm,stg, cdo,
      	        ctrl["NFOLD"],ctrl["vtmax"],ctrl["mcalc"],jlist,ctrl["S"],sd,σ)
      fvls[1:jsvd,sc] = tempa
      fvcs[1:jmsd,1:jsvd,sc] = tempe
      fqns[1:jsvd,:,sc] = tempq
   elseif ctrl["stages"] == 2
      #printstyled("hi wes\n",color=:cyan)
      tempa, tempe, tempt, tempq = twostg_calc(ctrl,prm,stg,cdo,
         ctrl["NFOLD"],ctrl["vtmax"],ctrl["mcalc"],ctrl["mmax"],jlist,ctrl["S"],sd,σ)
      fvls[1:jsvd,sc] = tempa
      fvcs[1:jmsd,1:jsvd,sc] = tempe
      fqns[1:jsvd,:,sc] = tempq
      tvcs[:,:,sc] = tempt
   else
      @warn "I can't diagonalize in this number of stages"
   end#if
   end#for
   println("yay energy levels are calculated!")
   #write energies to file
   if occursin("E",ctrl["RUNmode"])
      engwriter(molnam,ctrl,fvls,fqns)
   end
   scls = vcat(ser,cde)
   return fvls, fvcs, fqns, μs, prm, scls, stg, cdo, tvcs
end
function westersim(molnam::String,prm,ctrl,fvls,fvcs,fqns,μs,prms,scls,stg,ops,pcov,tvcs)
   #calculate transitions
#   if occursin("S",ctrl["RUNmode"])
   if pcov == nothing
      pcov = approxcovar(prm,permdeterm(scls,stg))
   end
   σcnt = σcount(ctrl["NFOLD"])
   jmax = ctrl["Jmax"]
   #@show pcov
   perm = permdeterm(scls,stg)
   #@show perm
   kbT = ctrl["TK"]*20836.61912 #MHz/K
   Qrt = sum(exp.(fvls ./ -kbT))/3
   finfrq = zeros(0,4)
   finqns = zeros(Int,0,12)
   #=
   if (ctrl["NFOLD"] > 0)&&(isodd(ctrl["NFOLD"]))
      σlst = [collect(1:σcnt) collect(1:σcnt)] .- 1
   elseif (ctrl["NFOLD"] > 0)&&(iseven(ctrl["NFOLD"]))
      σlst = [collect(1:σcnt) collect(1:σcnt); 1 σcnt; σcnt 1] .- 1
   else
      σlst = [0 0]
   end
   =#
   for sc in 1:σcnt
      σ = sc - 1
      vals = fvls[:,sc]
      vecs = fvcs[:,:,sc]
      quns = fqns[:,:,sc]
      if ctrl["stages"]==1
      J = derivcalc_all(ops,ctrl,perm,vecs,prm,scls,stg,σ)
      uncs = traerrs(J,pcov)
      else
         println("Sorry, the unc calculator isn't implemented"*
            " for this stage value ($(ctrl["stages"]))")
         uncs = zeros(size(vals))
      end
      #if σ==0
      #   @show (uncs)
      #end
      if ctrl["stages"]==1
      fr,qn = tracalc_nocat(μs,kbT,Qrt,ctrl,jmax,vals,vecs,quns,σ,
                           vals,vecs,quns,σ,uncs)
      elseif ctrl["stages"]==2
      fr,qn = tracalc_twostg(μs,kbT,Qrt,ctrl,jmax,vals,vecs,quns,σ,
                           vals,vecs,quns,σ,uncs,tvcs)
      else
         @warn "I can't diagonalize in this number of stages"
      end#if
      finfrq = vcat(finfrq,fr)
      finqns = vcat(finqns,qn)
   end
   #write transitions to file
   TraWriter(molnam, ctrl["S"], finfrq, finqns)
   return finfrq, finqns
#   end
end

function westerfit(molnam::String,ctrl::Dict{String,Any})
"""
   The fitter!
"""
   println("westerfit!")
   prm, ser = secordinp(molnam,String(ctrl["Irrep"]))
   μs, cdf, cdn, cde, cdo, stg = opinp(molnam,ctrl["NFOLD"])
   prm = vcat(prm,cdf)
   err = vcat(ser,cde)
   if occursin("F",ctrl["RUNmode"]) #Normal Fit behavior, overpowers check
      lines = readdlm("$molnam.lne", ',', Float64,comments=true,comment_char='#')
      linelength = (size(lines,1))
   else # Self-consistency check
      lines = readdlm("$molnam.cat", ',', Float64,comments=true,comment_char='#')
      lines = pred2lne(lines,ctrl["S"])
   end

   #determine the states
   linds, ofreqs, luncs = lineprep(lines,ctrl["NFOLD"],ctrl["S"],ctrl["vtmax"])
   #@show linds
   jlist = jlister(linds)
   #@show jlist
   #opt
#   println("Beginning optimization")
   outputinit(molnam,prm,err,linelength,ctrl)
if ctrl["stages"] == 1
   #tsrp, pcov, omcs, cfrqs, vals = lbmq_approx(ctrl,jlist,ofreqs,luncs,linds,prm,err,cdo,stg,molnam)
   tsrp, pcov, omcs, cfrqs, vals = lbmq_opttr(ctrl,jlist,ofreqs,luncs,linds,prm,err,cdo,stg,molnam)
elseif ctrl["stages"] == 2
   tsrp, pcov, omcs, cfrqs, vals = lbmq_2stg(ctrl,jlist,ofreqs,luncs,linds,prm,err,cdo,stg,molnam)
else
   @warn "I can't diagonalize in zero stages"
end#if
   #println(tsrp)
   #println(puncs)
   reswritter(molnam,lines,omcs,cfrqs)
   #println("New Parameter Vector:")
   #println("New Energy levels")
   #for n in 1:Nmax
   #   vals, vecs = rotdiag(Nmax,n,rotparams)
   #   println(vals)
   #end
   return tsrp, pcov
end

function westerfit(molnam::String)
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


#@time westerfit(ARGS[1])
