
using Dates
using DelimitedFiles
using FunctionWrappers
import FunctionWrappers: FunctionWrapper
using LinearAlgebra
import Printf: @sprintf
#import speed of light, plank, atomic mass, electron mass, fine structure constant
import PhysicalConstants.CODATA2022: c_0, k_B, h
using SparseArrays
using TOML
#@static if Sys.iswindows()
#@static if true
using WignerSymbols
wig3j(a,b,c,d,e,f) = wigner3j(Float64,a,b,c,d,e,f)
wig6j(a,b,c,d,e,f) = wigner6j(Float64,a,b,c,d,e,f)
#else
#   using WIGXJPFjl
#end

const NUMTYPE = Float64

#using JET
#using BenchmarkTools
#using ProfileView
include("@__DIR__/../psi.jl")
include("@__DIR__/../type.jl")
include("@__DIR__/../file_in.jl")
include("@__DIR__/../common.jl")
include("@__DIR__/../file_out.jl")
include("@__DIR__/../hamil.jl")
include("@__DIR__/../assign.jl")
include("@__DIR__/../ntop.jl")
include("@__DIR__/../transitions.jl")
include("@__DIR__/../lsq/opt-com.jl")
include("@__DIR__/../lsq/optimizer.jl")
include("@__DIR__/../derivs.jl")

include("@__DIR__/../ops/baseops.jl")
include("@__DIR__/../ops/dipoles.jl")

const csl::Float64 = (c_0 * 1e-4).val # MHz / cm⁻¹
const kb::Float64 = (k_B/h * 1e-6).val # MHz / K

BLAS.set_num_threads(Int(0.5*Sys.CPU_THREADS))
#@show Threads.nthreads()
#@warn "FUCK FUCK FUCK Ψ LENGTH IS MESSED UP. IT NEEDS TO KNOW ABOUT STAGES. FUCK ONE STAGE"
#@warn "There is a bug where vtmax & Jmax cause the fitter to break.
##This happens if they are too big/small relative to the line list"

if NUMTYPE <: Complex
   @warn "You have engaged C₁ mode. God have mercy on your soul & your runtimes"
end

function westereng(molnam::String, ctrl::Controls,ℋ::Vector{Term})::Eigs 
   wvs = Eigs(ctrl)
   jsσs = jσlister_full(ctrl.S,ctrl.Jmax, σcount(ctrl.NFOLD))
   H_calc(ctrl,wvs,ℋ,jsσs)
   if occursin("E",ctrl.RUNmode)
      engwriter(molnam, ctrl.Jmax, ctrl.S, ctrl.vtmax, wvs.rst.vals)
      println("yay! energy levels writen to $molnam.eng")
   end
   return wvs
end

function westersim(molnam::String, ctrl::Controls, μs::Vector{MuOp}, wvs::Eigs)
   σs = σcount(ctrl.NFOLD)
   frqs, inds = tracalc(ctrl,μs,wvs)
   writefreqs(molnam,ctrl,frqs,inds)
   return inds, frqs
end
function westerfit(molnam::String, ctrl::Controls, ℋ::Vector{Term})
   lins, qns = linereader(ctrl, molnam)
   ℋ, omc, cfrqs = opt_calc(molnam, ctrl, ℋ, lins)
   reswritter(molnam, qns, lins, omc, cfrqs)
   return ℋ
end

function westermain()
   println("Sorry about the name...")
   molnam = "test_sr"
   @time info, ctrl, ℋ, μs = inp_reader(molnam)
   if occursin("F", ctrl.RUNmode)
      println("westerfit!")
      @time wvs = westerfit(molnam, ctrl, ℋ)
   end # F
   if occursin("E", ctrl.RUNmode)||occursin("S", ctrl.RUNmode)
      @time wvs = westereng(molnam, ctrl, ℋ)
      if occursin("S", ctrl.RUNmode)
         @time westersim(molnam, ctrl, μs, wvs)
      end # S
   end # E / S
   return wvs
end


#Base.@ccallable function main()::Cint
function main()
   westermain(ARGS[1])
   return 0
end
