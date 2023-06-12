
using DelimitedFiles
using LinearAlgebra
using LinearAlgebra.BLAS
using Printf
using SparseArrays
using StaticArrays
using Base.Threads
using Dates
include("@__DIR__/../WIGXJPF.jl")
using .WIGXJPF
include("@__DIR__/../assign.jl")
include("@__DIR__/../common.jl")
include("@__DIR__/../filehandling.jl")
include("@__DIR__/../jacobi.jl")
include("@__DIR__/../optimizer.jl")
include("@__DIR__/../transitions.jl")
include("@__DIR__/../userdefop.jl")

BLAS.set_num_threads(Threads.nthreads())

const csl = 29979.2458

function westersim(molnam::String, prm, ctrl)
   println("westersim!")
   molnam = String(split(molnam,'.')[1])
   #molnam = replace(molnam, r".inp"=>"")
   #read input file
   sof, ser = secordinp(molnam)
   μs, cdf, cdn, cde, cdo, stg = opinp(molnam)
   if prm==nothing
      prm = vcat(sof,cdf)
   else
      prm[1:15] = sod2prep(prm[1:15])
   end
   #calculate energy levels
   σcnt = σcount(ctrl["NFOLD"])
   sd = Int(2.0*ctrl["S"]+1.0)
   jmin = 0.5*iseven(sd)
   jmax = ctrl["Jmax"]
   jlist = collect(Float64,jmin:jmax)
   vtd = Int(ctrl["vtmax"]+1)
   mcd = Int(2*ctrl["mcalc"]+(iseven(ctrl["NFOLD"]))+1)
   jfd = sd*Int(sum(2.0 .* jlist .+ 1.0))
   fvls = zeros(Float64,jfd*vtd,σcnt)
   fqns = zeros(Int,jfd*vtd,6,σcnt)
   fvcs = zeros(Float64,Int(sd*(2*jmax+2)*mcd),jfd*vtd,σcnt)
   @time @simd for sc in 1:σcnt
      σ = sc - 1
	   mcd = Int(2*ctrl["mcalc"]+(σtype(ctrl["NFOLD"],σ)==2)+1)
	   jmsd = Int(mcd*sd*(2*jmax+1))
	   jsvd = Int(jfd*vtd)
#       fvls[1:jvsd,σ], fvcs[1:jmsd,1:jsvd,σ], fqns[:,1:jvsd,σ] = tsrcalc(sof,cdf,cdo,
      tempa, tempe, tempq = tsrcalc(ctrl,prm,stg, cdo,
      	        ctrl["NFOLD"],ctrl["vtmax"],ctrl["mcalc"],jlist,ctrl["S"],sd,σ)
      fvls[1:jsvd,sc] = tempa
      fvcs[1:jmsd,1:jsvd,sc] = tempe
      fqns[1:jsvd,:,sc] = tempq
   end
   println("yay energy levels are calculated!")
   #write energies to file
   if occursin("E",ctrl["RUNmode"])
      engwriter(molnam,ctrl,fvls,fqns)
   end
   #calculate transitions
   if occursin("S",ctrl["RUNmode"])
      kbT = ctrl["TK"]*20836.61912 #MHz/K
      Qrt = sum(exp.(fvls ./ -kbT))/3
      finfrq = zeros(0,3)
      finqns = zeros(Int,0,12)
      @time for sc in 1:σcnt
         σ = sc - 1
         vals = fvls[:,sc]
         vecs = fvcs[:,:,sc]
         quns = fqns[:,:,sc]
         fr,qn = tracalc_nocat(μs,kbT,Qrt,ctrl,jmax,vals,vecs,quns,σ,vals,vecs,quns,σ)
         finfrq = vcat(finfrq,fr)
         finqns = vcat(finqns,qn)
      end
   #write transitions to file
   TraWriter(molnam, ctrl["S"], finfrq, finqns)
   end
end

function westerfit(molnam::String,ctrl::Dict{String,Any})
"""
   The fitter!
"""
   println("westerfit!")
   prm, ser = secordinp(molnam)
   μs, cdf, cdn, cde, cdo, stg = opinp(molnam)
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
   linds, ofreqs, luncs = lineprep(lines,ctrl["NFOLD"],ctrl["S"],0)
   jlist = jlister(linds)
   #println(linds)
   #opt
#   println("Beginning optimization")
   outputinit(molnam,prm,err,linelength)

   tsrp, puncs, omcs, cfrqs, vals = lbmq_opttr(ctrl,jlist,ofreqs,luncs,linds,prm,err,cdo,stg,molnam)
   #println(tsrp)
   #println(puncs)
   reswritter(molnam,lines,omcs,cfrqs)
   #println("New Parameter Vector:")
   #println("New Energy levels")
   #for n in 1:Nmax
   #   vals, vecs = rotdiag(Nmax,n,rotparams)
   #   println(vals)
   #end
   #write output file
   return tsrp
end

function westerfit(molnam::String)
   molnam = String(split(molnam,'.')[1])
   #read input file
   ctrl = ctrlinp(molnam)
   if occursin("T",ctrl["RUNmode"])
      westersim(molnam,nothing, ctrl)
      westerfit(molnam, ctrl)
   else
      if occursin("F",ctrl["RUNmode"])
         prm = westerfit(molnam, ctrl)
      else
         prm = nothing
      end
      if occursin("E", ctrl["RUNmode"])||occursin("S", ctrl["RUNmode"])
         westersim(molnam, prm, ctrl)
      end
   end
end


#@time westerfit(ARGS[1])
