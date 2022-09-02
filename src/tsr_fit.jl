"""
Welcome to westersim_fr.jl

this simplified version only uses 2nd order parameters for ease of testing

This program is a new, more group-theoretically driven version of westersim.
The fr stands for free rotor as this version does not use the belgi 2-stage
   approach and instead does everything in a single stage. The matrix is twice
   Wang transformed to improve eigenvector definition

interal parameters:
      1  2  3   4   5    6    7   8  9 10 11 12
prs =[A; B; C; Dab; Fr; Frρ; V3; ao; a; b; d; η]

input parameters:
      1  2  3  4  5   6   7    8    9    10  11
prs =[A; B; C; δ; F; V3; ϵzz; ϵxx; ϵyy; ϵxzb; η]

"""

using DelimitedFiles
#using Optim
using Printf
using LinearAlgebra
using LinearAlgebra.BLAS
using LinearAlgebra.LAPACK
#using LineSearches
using SparseArrays
include("/home/wes/files/westerfit/src/assign.jl")
include("/home/wes/files/westerfit/src/common.jl")
include("/home/wes/files/westerfit/src/filehandling.jl")
include("/home/wes/files/westerfit/src/hamiltonian.jl")
include("/home/wes/files/westerfit/src/intensities.jl")
include("/home/wes/files/westerfit/src/jacobi.jl")
#using .jacobi
include("/home/wes/files/westerfit/src/optimizer.jl")
#using Threads
include("/home/wes/files/westerfit/src/WIGXJPF.jl")
using .WIGXJPF

#using WignerSymbols
#wig3j(j1,j2,j3,m1,m2,m3) = Float64(wigner3j(j1,j2,j3,m1,m2,m3))
#wig6j(j1,j2,j3,j4,j5,j6) = Float64(wigner6j(j1,j2,j3,j4,j5,j6))

BLAS.set_num_threads(12)

apology = true

println("westerfit!")
print("Enter molecule name: \n")
#molnam = readline()
#molnam = "mepho-m"
#include("$molnam.par")
#println(molnam)
#println("Nmax=$Nmax")
if apology==true
   println("Sorry about the name...")
end

################################################################################
############                       westersim                        ############
################################################################################
function tsrdiag(pr,j,s,nmax,mcalc,mmax,σ)
   U = ur(j,s,mcalc)
   if σ==zero(σ)
      U = ur(j,s,mcalc)*ut(mcalc,j,s)
      H = Matrix(U*Htsrmat2(pr,j,s,mcalc,σ)*U)
      #H = Matrix(Htsrmat2(pr,j,s,mcalc,σ))
   else
      H = Matrix(Htsrmat2(pr,j,s,mcalc,σ))
   end
   H, rvec = jacobisweep2(H,Int(floor((j+s+6*mcalc)/8)))
   vals, vecs = LAPACK.syev!('V', 'U', H)
   qns, vals, vecs = assign(j,s,σ,mcalc,vals,vecs,rvec)
   return qns, vals, vecs
end

function tsrcalc(prm,s,nmax,mcalc,mmax,σ)
   #prm = PAM2RAM(prm)
   #prm = paramprep(prm)
   #prm = paramshift(inprms)
   jmin =  0.5*iseven(Int64(2*s+1))
   jmax = nmax - s
   jarray = collect(Float64,jmin:jmax)
   jd = Int((2*s+1)*sum(2 .* jarray .+ 1))
   outvals = zeros(Float64,Int(jd*(2*mmax+1)))
   outvecs = zeros(Float64,Int((2*s+1)*(2*jmax+2)*(2*mcalc+1)),Int(jd*(2*mmax+1)))
   outqns = zeros(Float64,Int(jd*(2*mmax+1)),6)
   Threads.@threads for i in 1:length(jarray)
      j = jarray[i]
      sind, find = jinds(j,s,mmax)
      tqns, tvals, tvecs = tsrdiag(prm,j,s,nmax,mcalc,mmax,σ)
      si = findfirst(isequal(-mmax*NFOLD+σ),tqns[:,5])
      fi = findlast(isequal(mmax*NFOLD+σ),tqns[:,5])
      outvals[sind:find] = tvals[si:fi]
      outvecs[1:Int((2*s+1)*(2*j+1)*(2*mcalc+1)),sind:find] = tvecs[:,si:fi]
      outqns[sind:find,:] = tqns[si:fi,:]
   end
   return outqns, outvals, outvecs
end

function westersim(prm,s,nmax,mcalc,mmax)
   σs = σcount(NFOLD)
   jmin =  0.5*iseven(Int64(2*s+1))
   jmax = nmax - s
   jarray = collect(Float64,jmin:jmax)
   jd = Int((2*s+1)*sum(2 .* jarray .+ 1))
   vals = zeros(Float64,Int(jd*(2*mmax+1)),σs)
   vecs = zeros(Float64,Int((2*s+1)*(2*jmax+2)*(2*mcalc+1)),Int(jd*(2*mmax+1)),σs)
   qns = zeros(Float64,Int(jd*(2*mmax+1)),6,σs)
   @time qns[:,:,1], vals[:,1], vecs[:,:,1] = tsrcalc(parameters,s,nmax,mcalc,mmax,0)
   @time transitions = tracalc(nmax,s,mcalc,0,qns[:,:,1],vals[:,1],vecs[:,:,1])
   for σ in 1:(σs-1)
      @time qns[:,:,σ+1], vals[:,σ+1], vecs[:,:,σ+1] = tsrcalc(parameters,s,nmax,mcalc,mmax,σ)
      @time temptrns = tracalc(nmax,s,mcalc,σ,qns[:,:,σ+1],vals[:,σ+1],vecs[:,:,σ+1])
      transitions = vcat(transitions,temptrns)
   end
   return transitions, qns, vals, vecs
end

################################################################################
############                  westerfit                  ############
################################################################################
function xtxsolve(A)
   x = zeros(size(A))
   a,b = eigen(A)
   x = √(diagm(a))*transpose(b)
   return x
end

function jlister(inds)
   #finds all the unique J & σ pairs
   js = vcat(inds[:,1],inds[:,4])
   σs = vcat(inds[:,2],inds[:,5])
   temp = fill((0,0),size(js))
   for i in 1:size(js)[1]
      temp[i] = (js[i],σs[i])
   end
   temp = unique(temp)
   jsσs = zeros(Int,size(temp)[1],2)
   jsσs[:,1] = (x->x[1]).(temp)
   jsσs[:,2] = (x->x[2]).(temp)
   jsσs = jsσs[sortperm(jsσs[:,1]),:]
#   jsσs = jsσs[sortperm(jsσs[:,2])]
   return jsσs
end

function limeigcalc(jlist,inds,rp)
"""
This is the 'limited eigen calculator'. It operates similarly to the tsrdiag
   function except for only the J-σ pairs from the input line list rather than
   covering every pairing.
"""
   #rp = paramprep(rp)
   jmax = 0.5*maximum(jlist[:,1])
   #println(jlist)
   σs = σcount(NFOLD)
   jd = (2*S+1)*sum(2 .* collect((0.5*isodd(2*S)):jmax) .+ 1)
   #evals = zeros(Float64,Int(jd*(2*mmax+1)),σs)
   #evecs = zeros(Float64,Int((2*S+1)*(2*jmax+1)*(2*mcalc+1)),Int(jd*(2*mmax+1)),σs)
   #eqns = zeros(Float64,Int(jd*(2*mmax+1)),6,σs)
   evals = zeros(Float64,Int(jd*(2*mmax+1)),σs)
   evecs = zeros(Float64,Int((2*S+1)*(2*jmax+2)*(2*mcalc+1)),Int(jd*(2*mmax+1)),σs)
   eqns = zeros(Float64,Int(jd*(2*mmax+1)),6,σs)
   Threads.@threads for i in 1:size(jlist)[1]
      j = 0.5*jlist[i,1]
      σ = jlist[i,2]
      sind, find = jinds(j,S,mmax)
      tqns, tvals, tvecs = tsrdiag(rp,j,S,nmax,mcalc,mmax,σ)
      si = findfirst(isequal(-mmax*NFOLD+σ),tqns[:,5])
      fi = findlast(isequal(mmax*NFOLD+σ),tqns[:,5])
      evals[sind:find,σ+1] = tvals[si:fi]
      evecs[1:Int((2*S+1)*(2*j+1)*(2*mcalc+1)),sind:find,σ+1] = tvecs[:,si:fi]
      eqns[sind:find,:,σ+1] = tqns[si:fi,:]
   end
   return evals, evecs#, eqns
end
#take difference and determine rms



function dogleg(βlm,βsd,t,Δ)
"""
This enforces a trust region on our steps
"""
   βf = zeros(size(βlm))
   nβlm = norm(βlm)
   nβsd = norm(βsd)
   if nβlm ≤ Δ #lm step inside trust region
      βf = βlm
   elseif (nβlm>Δ)&&(nβsd≤Δ) #only SD step inside trust region
      βlts = βlm - t.*βsd
      s = (Δ-abs(t)*nβsd)/norm(βlts) #this is an aproximate solution ≤ the actual
      βf = t*βsd + s.*(βlts)
   else#if (norm(βlm)>Δ)&&(norm(βsd)>Δ) #both outside trust region
      βf = (Δ/nβsd) .* βsd
   end
   return βf
end

function westerfit_handcoded()
   global mmax=0
   global mcalc=1
   global S=0.5
#   tsrparams = PAM2RAM(parameters)
   tsrparams = parameters
   #println(tsrparams)
   #read the lines
#   lines = readdlm("$molnam.lne", ',', Float64)
   lines = testlines
   #determine the states
   linds, ofreqs, uncs = lineprep(lines,S,mmax)
   println(linds)
   jlist = jlister(linds)
   global nmax= S + 0.5*maximum(jlist[:,1])
   #opt
   scales = trsscales
   λ = 1.0E+4
#   println("Beginning optimization")
   tsrp, vals = lbmq_opttr(jlist,ofreqs,uncs,linds,tsrparams,scales,λ)
   #println("New Parameter Vector:")
   #println(tsrp)
   println(tsrp - oparameters)
   #println("New Energy levels")
   #for n in 1:Nmax
   #   vals, vecs = rotdiag(Nmax,n,rotparams)
   #   println(vals)
   #end
   #write output file
end

################################################################################
############                    test parameters                     ############
################################################################################

const csl = 29979.2458 #MHz to cm-1 converter
F = 292574.52
V3 = 4195821.0
A = 82756.97
B = 10145.212
C = 9426.952
#δ = 0.45
ρ = 0.0513448477
Dab = -3716.8
ϵzz = 254.3618
ϵxx = 385.4886
ϵyy = -415.9114
ϵxz = 345.6789
ΔN = 0.01188
η = 0.0
μa = 1.0
μb = 0.5
μ = [μa 0.008; μb 0.006; 0.0 -0.034]

#parameters = [ A;   B;   C;   δ;   F;   V3; ϵzz; ϵxx; ϵyy; ϵxz;  η]
#parameters = [ A+F*ρ^2;   B;   C;   Dab;   F; -ρ*F;  V3; ϵzz; ϵxx; ϵyy; ϵxz;  η; ΔN]
#trsscales = ones(Float64,size(parameters))

#trsscales =  [0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01]
#trsscales = fill(0.1,13)

"""
testlines = [
 1.5  1  0  1  0  0  0.5  0  0  0  0  0  19419.0217  0.4; #2
 0.5  1  0  1  0  0  0.5  0  0  0  0  0  19450.6917  0.4; #1
 2.5  2  0  2  0  0  0.5  1  0  1  0  0  38819.1124  0.4; #
 2.5  2  0  2  0  0  1.5  1  0  1  0  0  38847.8178  0.4; #7
 1.5  2  0  2  0  0  0.5  1  0  1  0  0  38865.0721  0.4; #5
 1.5  2  0  2  0  0  1.5  1  0  1  0  0  38894.2759  0.4; #6
 2.5  3  0  3  0  0  1.5  2  0  2  0  0  58287.5536  0.4; #9
 3.5  3  0  3  0  0  2.5  2  0  2  0  0  58272.0583  0.4; #11
 1.5  1  0  1  1  1  0.5  0  0  0  1  1  19104.1433  0.4; #16
 0.5  1  0  1  1  1  0.5  0  0  0  1  1  19137.5883  0.4; #15
 2.5  2  0  2  1  1  1.5  1  0  1  1  1  38250.8564  0.4; #21
 1.5  2  0  2  1  1  0.5  1  0  1  1  1  38263.2012  0.4; #19
 1.5  2  0  2  1  1  1.5  1  0  1  1  1  38295.7924  0.4; #20
 2.5  2  0  2  1  1  0.5  1  0  1  1  1  38228.7521  0.4]
"""

INTTHRESHOLD = .00000
#=
parameters:
     1  2  3   4   5  6   7   8  9 10 11 12
prs =[A; B; C; Dab; Fr; ρ; V3; ao; a; b; d; η]

λs.vectors[abs.(λs.vectors) .< 1E-9] .= 0.0
vecs[abs.(vecs) .< 1E-9] .= 0.0
=#
function westerenergies(prm,s,nmax,mcalc,mmax)
   @time amq, ama, amv = tsrcalc(parameters,s,nmax,mcalc,mmax,0)
   EngWriter(ama, amq, mmax, 0)
   println("E states!")
   @time emq, ema, emv = tsrcalc(parameters,s,nmax,mcalc,mmax,1)
   EngWriter(ema, emq, mmax, 1)
   qns  = cat(amq,emq, dims=3)
   vals = cat(ama,ema, dims=3)
   vecs = cat(amv,emv, dims=3)
   return qns, vals, vecs
end


global NFOLD=3
global TK = 25.0
const KB = 2.083661912E+4 #MHz/K
mc = 1
nm = 3

#=
A      =      0.35150076793036794*29979.2458
B      =      0.15346996670417229*29979.2458
C      =      0.10623368397991501*29979.2458
Dab    =     -0.03670915325930996*29979.2458
V3     =    359.14177475235487691*29979.2458
ρ      =      0.06376911062071215
F      =      5.64133778441192124*29979.2458
parameters = [ A+F*ρ^2;   B;   C; 0*Dab;   F; ρ*F;  V3; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
=#

parameters = [  A;   B;   C; Dab;   F; ρ*F;  V3; ϵzz; ϵxx; ϵyy; ϵxz; 0.0;  ΔN]
trsscales = [ 1.0; 1.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0]


molnam = "hirota"
@time transitions, qns, vals, vecs = westersim(parameters,0.5,nm,mc,0)
testlines = (pred2lne(transitions))
println(size(testlines))
oparameters = copy(parameters)

pert = 0.001*(0.5 .- rand(Float64,size(parameters))).*trsscales.*parameters
println(pert)
#pert = trsscales
parameters .+= pert
#orgparams = parameters
@time westerfit_handcoded()
