
using DelimitedFiles
using LinearAlgebra
using LinearAlgebra.BLAS
using Printf
using SparseArrays
using StaticArrays
using Base.Threads
include("@__DIR__/../WIGXJPF.jl")
using .WIGXJPF
include("@__DIR__/../common.jl")
include("@__DIR__/../filehandling.jl")
include("@__DIR__/../jacobi.jl")
include("@__DIR__/../optimizer.jl")
include("@__DIR__/../transitions.jl")
include("@__DIR__/../userdefop.jl")

BLAS.set_num_threads(Threads.nthreads())

const csl = 29979.2458

function westersim(molnam::String, ctrl)
   println("westersim!")
   molnam = replace(molnam, r".inp"=>"")
   #read input file
   prm, ser = secordinp(molnam)
   μs, cdf, cdn, cde, cdo, stg = opinp(molnam)
   prm = vcat(prm,cdf)
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
      tempa, tempe, tempq = tsrcalc(prm,stg, cdo,
      	        ctrl["NFOLD"],ctrl["vtmax"],ctrl["mcalc"],jlist,ctrl["S"],sd,σ)
      fvls[1:jsvd,sc] = tempa
      fvcs[1:jmsd,1:jsvd,sc] = tempe
      fqns[1:jsvd,:,sc] = tempq
   end
   println("yay energy levels are calculated!")
   #write energies to file
   if occursin("E",crtl["RUNmode"])
      EngWriter(ctrl,fvls,fqns)
   end
   #calculate transitions
   if occursin("S",crtl["RUNmode"])
      kbT = ctrl["TK"]*20836.61912 #MHz/K
      Qrt = sum(exp.(fvls))
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

function westerfit(molnam::String)
   molnam = replace(molnam, r".inp"=>"")
   #read input file
   ctrl = ctrlinp(molnam)
   if occursin("E",crtl["RUNmode"])||occursin("S",crtl["RUNmode"])
      westersim(molnam, ctrl)
   end
   if occursin("F",ctrl["RUNmode"])
      westerfit(molnam, ctrl)
   end
end


#h = hsr(ones(12),0.5,0.5)
#println(h)

function westerfit(molnam::String,ctrl)
"""
   The fitter!
"""
   println("westerfit!")
   prm, ser = secordinp(molnam)
   μs, cdf, cdn, cde, cdo, stg = opinp(molnam)
   prm = vcat(prm,cdf)
   err = vcat(ser,cde)
   lines = readdlm("$molnam.lne", ',', Float64)
   #determine the states
   linds, ofreqs, uncs = lineprep(lines,ctrl["NFOLD"],ctrl["S"],0)
   #println(linds)
   jlist = jlister(linds)
   #global nmax = S + 0.5*maximum(jlist[:,1])
   #opt
#   println("Beginning optimization")
   tsrp, vals = lbmq_opttr(ctrl,jlist,ofreqs,uncs,linds,prm,err,cdo,stg)
   #println("New Parameter Vector:")
   println(tsrp)
   #println("New Energy levels")
   #for n in 1:Nmax
   #   vals, vecs = rotdiag(Nmax,n,rotparams)
   #   println(vals)
   #end
   #write output file
end

@time westerfit(ARGS[1])
