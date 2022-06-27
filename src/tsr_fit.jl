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
using Optim
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
#using Threads
include("/home/wes/files/westerfit/src/WIGXJPF.jl")
using .WIGXJPF

BLAS.set_num_threads(8)

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

global NFOLD=3
global TK = 25.0

################################################################################
############                       westersim                        ############
################################################################################


function tsrdiag(pr,j,s,nmax,mcalc,mmax,σ)
   U = ur(j,s,mcalc)
   if σ==0||σ==0.0
      U = ur(j,s,mcalc)*ut(mcalc,j,s)
      H = Matrix(U*Htsr(pr,j,s,mcalc,σ)*U)
   else
      #U = ur(j,s,mcalc)
      #H = Matrix(U*Htsr(pr,j,s,mcalc,σ)*U)
      H = Matrix(Htsr(pr,j,s,mcalc,σ))
   end
   H, rvec = jacobisweep(H,floor((j+s)/2 - 1 + mcalc*σ/3 ))
   vals, vecs = LAPACK.syev!('V', 'U', H)
   #if (Int(j)==4)&&(Int(σ)==1.0)
   #   println(vals ./ csl)
   #   println(vecs)
   #end
   qns, vals, vecs = assign(nmax,j,s,σ,mcalc,mmax,vals,vecs,rvec)
   return qns, vals, vecs
end

function tsrcalc(prm,s,nmax,mcalc,mmax,σ)
   #prm = PAM2RAM(prm)
   #prm = paramprep(prm)
   #prm = paramshift(inprms)
   if isodd(Int64(2*s+1))
      jmin = 0.0
   else
      jmin = 0.5
   end
   jmax = nmax - s
   jarray = collect(Float64,jmin:jmax)
   jd = Int((2*s+1)*sum(2 .* jarray .+ 1))
   outvals = zeros(Float64,Int(jd),2*mmax+1)
   outvecs = zeros(Float64,Int((2*s+1)*(2*jmax+2)*(2*mcalc+1)),Int(jd*(2*mmax+1)),2*mmax+1)
   outqns = zeros(Float64,jd,6,2*mmax+1)
   Threads.@threads for i in 1:length(jarray)
      j = jarray[i]
      sind, find = jinds(j,s)
      tqns, tvals, tvecs = tsrdiag(prm,j,s,nmax,mcalc,mmax,σ)
      si = 1
      f0 = Int((2*s+1)*(2*j+1))
      fi = f0
      for m in 1:(2*mmax+1)
         outvals[sind:find,m] = tvals[si:fi]
         outvecs[1:Int((2*s+1)*(2*j+1)*(2*mcalc+1)),sind:find,m] = tvecs[:,si:fi]
         outqns[sind:find,:,m] = tqns[si:fi,:]
         si += f0
         fi += f0
      end
   end
   mperm = sortperm(collect(-mmax:mmax),by=abs)
   outvals = outvals[:,mperm]
   outvecs = outvecs[:,:,mperm]
   outqns = outqns[:,:,mperm]
   return outqns, outvals, outvecs
end

function westersim(prm,s,nmax,mcalc,mmax)
   @time amq, ama, amv = tsrcalc(parameters,s,nm,mc,0,0)
   #EngWriter(ama, amq, 0)
   @time atransitions = tracalc(nm,s,mc,0,amq,ama,amv)
   println("E states!")
   @time emq, ema, emv = tsrcalc(parameters,s,nm,mc,0,1)
   @time etransitions = tracalc(nm,s,mc,1,emq,ema,emv)
   #(nmax,s,mcalc,σ,qns,vals,vecs)
   #EngWriter(ema, emq,1)
   qns  = cat(amq,emq, dims=3)
   vals = cat(ama,ema, dims=3)
   vecs = cat(amv,emv, dims=3)
   transitions = vcat(atransitions,etransitions)
   return transitions, qns, vals, vecs
end

################################################################################
############                  westerfit                  ############
################################################################################

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
   evals = zeros(Float64,Int(jd*(2*mmax+1)),σs)
   evecs = zeros(Float64,Int((2*S+1)*(2*jmax+1)*(2*mcalc+1)),Int(jd*(2*mmax+1)),σs)
   eqns = zeros(Float64,Int(jd*(2*mmax+1)),6,σs)
   Threads.@threads for i in 1:size(jlist)[1]
      j = 0.5*jlist[i,1]
      σ = jlist[i,2]
      sind, find = jinds(j,S)
      tqns, tvals, tvecs = tsrdiag(rp,j,S,nmax,mcalc,mmax,σ)
      evals[sind:find,σ+1] = tvals
      evecs[1:Int((2*S+1)*(2*j+1)*(2*mcalc+1)),sind:find,σ+1] = tvecs
      eqns[sind:find,:,σ+1] = tqns
   end
   return evals, evecs#, eqns
end
#take difference and determine rms
function cost(vals,inds,ofreqs,weights)
   cfreqs = zeros(size(ofreqs))
   for i in 1:size(cfreqs)[1]
      cfreqs[i] = vals[inds[i,3],inds[i,2]+1] - vals[inds[i,6],inds[i,5]+1]
   end
   omc = ofreqs - cfreqs
   cst = sum((omc ./weights) .^2)
   rms = sqrt(cst/length(omc))
   return cst
end
function rmscalc(vals,inds,ofreqs)
   cfreqs = zeros(size(ofreqs))
   for i in 1:size(cfreqs)[1]
      cfreqs[i] = vals[inds[i,3],inds[i,2]+1] - vals[inds[i,6],inds[i,5]+1]
   end
   omc = ofreqs - cfreqs
   rms = sqrt(sum(omc .^ 2)/length(omc))
   return rms, omc
end
#construct jacobian
function anaderiv(j,s,σ,vec,rp,rpid)
   rp = zeros(Float64, length(rp)+1)
   rp[rpid] = 1.0
   U = ur(j,s,mcalc)
   if σ==0
      U *= ut(mcalc,j,s)
   end
   mat = Matrix(U*Htsr(rp,j,s,mcalc,σ)*U)
   out = transpose(vec)*mat*vec
   return out[1]
end
function build_jcbn(inds,vecs,params)
"""
This builds the Jacobian based on the Hellmann–Feynman theorem.
"""
   #fitprm = params[perm]
   jcbn = zeros(Float64,size(inds)[1],length(params))
   Threads.@threads for a in 1:size(inds)[1]
      ju = 0.5*inds[a,1]
      jl = 0.5*inds[a,4]
      σu = inds[a,2]
      σl = inds[a,5]
      vecu = vecs[1:Int((2*S+1)*(2*ju+1)*(2*mcalc+1)),inds[a,3],σu+1]
      #println(size(vecu))
      vecl = vecs[1:Int((2*S+1)*(2*jl+1)*(2*mcalc+1)),inds[a,6],σl+1]
      #println(vecs[:,inds[a,3],σu+1])
      #println(transpose(vecl)*vecl)
      for b in 1:length(params)
         jcbn[a,b] = anaderiv(ju,S,σu,vecu,params,b) - anaderiv(jl,S,σl,vecl,params,b)
      end
   end
   return jcbn
end

function dHdA(j,s,σ,vec,rp)
   drp = zeros(Float64, length(rp)+1)
   A = rp[1]
   F = rp[5]
   #dAeff/dAp
   drp[1] = (F^2) / ((F - A)^2)
   #dFr/dAp
   drp[5] = (F^2) / ((F - A)^2)
   #dFρ/dAp
   drp[6] = (F^2) / ((F - A)^2)
   if σ==0
      U = ur(j,s,mcalc)
      U *= ut(mcalc,j,s)
   else
      mat = Matrix(Htsr(drp,j,s,mcalc,σ))
   end
   out = transpose(vec)*mat*vec
   return out[1]
end
function dHdF(j,s,σ,vec,rp)
   drp = zeros(Float64, length(rp)+1)
   A = rp[1]
   F = rp[5]
   #dAeff/dF
   drp[1] = -(A^2) / ((F - A)^2)
   #dFr/dF
   drp[5] = (F^2 - 2.0A*F) / ((F - A)^2)
   #dFρ/dF
   drp[6] = -(A^2) / ((F - A)^2)
   if σ==0
      U = ur(j,s,mcalc)
      U *= ut(mcalc,j,s)
   else
      mat = Matrix(Htsr(drp,j,s,mcalc,σ))
   end
   out = transpose(vec)*mat*vec
   return out[1]
end
function build_jcbn2(inds,vecs,params)
"""
This builds the Jacobian based on the Hellmann–Feynman theorem. This is a
   modified version that uses particular expressions for F and A due to their
   direct couplings through ρ
"""
   #fitprm = params[perm]
   jcbn = zeros(Float64,size(inds)[1],length(params))
   for a in 1:size(inds)[1]
      ju = 0.5*inds[a,1]
      jl = 0.5*inds[a,4]
      σu = inds[a,2]
      σl = inds[a,5]
      vecu = vecs[1:Int((2*S+1)*(2*ju+1)*(2*mcalc+1)),inds[a,3],σu+1]
      vecl = vecs[1:Int((2*S+1)*(2*jl+1)*(2*mcalc+1)),inds[a,6],σl+1]
      jcbn[1,b] = dHdA(ju,S,σu,vecu,params,b) - dHdA(jl,S,σl,vecl,params,b)
      jcbn[2,b] = anaderiv(ju,S,σu,vecu,params,b) - anaderiv(jl,S,σl,vecl,params,b)
      jcbn[3,b] = anaderiv(ju,S,σu,vecu,params,b) - anaderiv(jl,S,σl,vecl,params,b)
      jcbn[4,b] = anaderiv(ju,S,σu,vecu,params,b) - anaderiv(jl,S,σl,vecl,params,b)
      jcbn[5,b] = dHdF(ju,S,σu,vecu,params,b) - dHdF(jl,S,σl,vecl,params,b)
      for b in 6:length(params)
         jcbn[a,b] = anaderiv(ju,S,σu,vecu,params,b) - anaderiv(jl,S,σl,vecl,params,b)
      end
   end
   return jcbn
end

function lbmq_step(jcbn, weights, omc, lbmq_prm)
"""
This should be the Levenberg-Marquadt step. This solves (JᵗWJ+λI)Δβ = (JᵗW)Δy
   for Δβ. Where J is the Jacobian, W is the weights, λ is the
   Levenberg-Marquadt parameter, and Δy is the omcs. This returns the step, Δβ.
"""
   jtw = transpose(jcbn)*weights
   β = zeros(size(jcbn)[2])
   jtj = jtw*jcbn
   A = Hermitian(jtj + lbmq_prm*Diagonal(jtj))
   A = cholesky(A)
   X = -jtw*omc
   β = ldiv!(β, A, X)
   β = β #./ norm(β.^2)
   return β
end
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
function acc_lbmq_step(jcbn, weights, omc, lbmq_prm)
"""
This should be the Geodesic Accelerated Levenberg-Marquadt step.
"""
   g = transpose(jcbn)*jcb + UniformScaling(lbmq_prm)
   ∇C = transpose(jcbn)*omc
   v = -inv(g)*∇C
   return Δx
end

function solvcalc(rprms,jlist,linds,ofreqs)
"""
This is the objective function for interfacing with other packages like Optim.jl
   It takes the parameters, the j σ pairs, line indices, and observed freqs.
   It then returns the rms. It can structured as x ->
   solvcalc(x,jlist,linds,ofreqs) for nice interfacing.
"""
   vals,vecs = limeigcalc(jlist, linds, rprms)
   rms, omc = rmscalc(vals, linds, ofreqs)
   println(rms)
   return rms
end
function solvcalc0(rprms,jlist,linds,ofreqs)
"""
This is the objective function for interfacing with other packages like Optim.jl
   It takes the parameters, the j σ pairs, line indices, and observed freqs.
   It then returns the rms. It can structured as x ->
   solvcalc(x,jlist,linds,ofreqs) for nice interfacing. It also prints the
   line indices, eigenvalues, and omcs for troubleshooting
"""
   vals,vecs = limeigcalc(jlist, linds, rprms)
   rms, omc = rmscalc(vals, linds, ofreqs)
   println(linds)
   println(vals)
   println(omc)
   return rms
end
function jcbn!(G,rps)
   G = build_jcbn(linds,vecs,rps)
   return G
end
function lbmq_opt(nlist,ofreqs,uncs,inds,params,scales,λ)
   vals,vecs = limeigcalc(nlist, inds, params)
   rms, omc = rmscalc(vals, inds, ofreqs)
   println("Initial RMS = $rms")
   counter = 0
   goal = sum(uncs)/length(uncs)
   weights = diagm(0=>(uncs .^ -1))
   THRESHOLD = 1.0E-06
   LIMIT = 100
   while true
      #println("Building Jacobians")
      jcbn = build_jcbn(inds,vecs,params)
      #println("Performing Levenberg-Marquadt Step")
      adjst = lbmq_step(jcbn,weights,omc,λ)
      adjst .*= scales
      #back up parameters
      newparams = params .+ adjst
      #recalculate RMS
      vals,vecs = limeigcalc(nlist, inds, newparams)
      nrms, omc = rmscalc(vals, inds, ofreqs)
      if nrms ≤ rms
         #accept step and decrease λ
         λ = λ/3.0
         params = newparams
         rms = nrms
      else #nrms > rms
         #reject step due to RMS increase
         λ = λ*2.0
         params = newparams
         rms = nrms
#      else
#         λ = λ
      end
      counter += 1
      check = abs(nrms-rms)/rms
      #println(check)
      println("After $counter interations, RMS = $rms, λ = $λ")
      if (rms ≤ goal)#||(check < THRESHOLD)||#&&(counter > 1)
         println("A miracle has come to pass. The fit has converged")
         break
      elseif counter ≥ LIMIT
         println("Alas, the iteration count has exceeded the limit")
         break
      else
         #write update to file
      end
   end
   return params, vals
end
function westerfit_handcoded()
#   tsrparams = PAM2RAM(parameters)
   tsrparams = parameters
   #read the lines
#   lines = readdlm("$molnam.lne", ',', Float64)
   lines = testlines
   #determine the states
   linds, ofreqs, uncs = lineprep(lines)
   jlist = jlister(linds)
   global mmax=0
   global mcalc=2
   global S=0.5
   global nmax= S + 0.5*maximum(jlist[:,1])
   #opt
   scales = trsscales
   λ = 0.1E-50
#   println("Beginning optimization")
   tsrp, vals = lbmq_opt(jlist,ofreqs,uncs,linds,tsrparams,scales,λ)
   println("New Parameter Vector:")
   println(tsrp)
   println("New Energy levels")
   #for n in 1:Nmax
   #   vals, vecs = rotdiag(Nmax,n,rotparams)
   #   println(vals)
   #end
   #write output file
end

function westerfit_optim()
"""
This version of westerfit uses the Optim package currently with the Accelerated
   Gradient Descent. It is profoundly slow but does seem to approach convergence.
   Best to allow it to run bit by bit. It is programed to only run for about
   600 seconds at a time. It will then print the new parameters so it can be
   restarted using those values.
"""
   x0 = parameters
   #read the lines
#   lines = readdlm("$molnam.lne", ',', Float64)
   lines = testlines
   #determine the states
   linds, ofreqs, uncs = lineprep(lines)
   #println(linds)
   goal = sum(lines[:,end])/length(lines[:,end])
   jlist = jlister(linds)
   global jlist=jlist
   global linds = linds
   global ofreqs = ofreqs
   global mmax=0
   global mcalc=8
   global S=0.5
   global nmax= S + 0.5*maximum(jlist[:,1])
#   scales = vcat(ones(3),zeros(1))
#   println("Beginning optimization")
   init_rms = solvcalc0(x0,jlist,linds,ofreqs)
   println("Initial RMS = $init_rms")
   lower = [82755.0, -1.0E-37, -1.0E-37, -1.0E+6, 1.0E-37, -1.0E+09, -1.0E+09,
             254.3617999, -1.0E+09, -1.0E+09, -1.0E-37, -1.0E-37]
   upper = [82758.0, 1.0E+09, 1.0E+09, 1.0E+6, 1.0E+12, 1.0E+09, 1.0E+09, 254.3618001,
             1.0E+09, 1.0E+09, 1.0E-37, 1.0E-37]
#   res = optimize(solvcalc, tsrparams, ParticleSwarm(; n_particles=12))
#   inner_optimizer = GradientDescent(linesearch=LineSearches.BackTracking(order=3))
#   res = optimize(solvcalc, lower, upper, tsrparams, Fminbox(inner_optimizer))
#         Optim.Options(iterations=100))
   res = Optim.optimize(x->solvcalc(x,jlist,linds,ofreqs), x0, AcceleratedGradientDescent(),
   Optim.Options(time_limit=600))#,
   println(res)
   println("Final RMS = ", Optim.minimum(res), " MHz")
   rotparams = Optim.minimizer(res)
   println("New Parameter Vector:")
   println(rotparams)
#   println("New Energy levels")
#   for j in 1:nmax
#      qns, vals, vecs = tsrdiag(tsrparams,j,s,nmax,mcalc,0,σ)
#      println(vals)
#   end
   #write output file
end

################################################################################
############               test parameters                ############
################################################################################



const csl = 29979.2458 #MHz to cm-1 converter
F = 292574.52
V3 = 4195821.0
A = 82756.97
B = 10145.212
C = 9426.952
#δ = 0.45
ρ = 0.513448477
Dab = -3716.8
ϵzz = 254.3618
ϵxx = 385.4886
ϵyy = -415.9114
ϵxz = 0.0
ΔN = 0.01188
η = 0.0
μa = 1.0
μb = 0.5

#parameters = [ A;   B;   C;   δ;   F;   V3; ϵzz; ϵxx; ϵyy; ϵxz;  η]
parameters = [ A-F*ρ^2;   B;   C;   Dab;   F; ρ*F;  V3; ϵzz; ϵxx; ϵyy; ϵxz;  η; ΔN]
#trsscales =  [0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01]
#parameters = [83316.10066771678, 14497.783454870654, 4733.339691158756,
#      -1.1414038835204492, 292650.6595414094, 0.0, 4.195803378991427e6,
#      306.43192861068763, 297.11625650644,
#      -22.797297784556505, 0.0, 0.0, 33.90230680364543]

trsscales = fill(1.0,12)

"""
testlines = [
 1.5  1  0  1  0  0  0.5  0  0  0  0  0  19419.0217  0.4; 2
 0.5  1  0  1  0  0  0.5  0  0  0  0  0  19450.6917  0.4; 1
 2.5  2  0  2  0  0  0.5  1  0  1  0  0  38819.1124  0.4;
 2.5  2  0  2  0  0  1.5  1  0  1  0  0  38847.8178  0.4; 7
 1.5  2  0  2  0  0  0.5  1  0  1  0  0  38865.0721  0.4; 5
 1.5  2  0  2  0  0  1.5  1  0  1  0  0  38894.2759  0.4; 6
 2.5  3  0  3  0  0  1.5  2  0  2  0  0  58287.5536  0.4; 9
 3.5  3  0  3  0  0  2.5  2  0  2  0  0  58272.0583  0.4; 11
 1.5  1  0  1  1  1  0.5  0  0  0  1  1  19104.1433  0.4; 16
 0.5  1  0  1  1  1  0.5  0  0  0  1  1  19137.5883  0.4; 15
 2.5  2  0  2  1  1  1.5  1  0  1  1  1  38250.8564  0.4; 21
 1.5  2  0  2  1  1  0.5  1  0  1  1  1  38263.2012  0.4; 19
 1.5  2  0  2  1  1  1.5  1  0  1  1  1  38295.7924  0.4; 20
 2.5  2  0  2  1  1  0.5  1  0  1  1  1  38228.7521  0.4]
"""

INTTHRESHOLD = .000000
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

mc = 7
nm = 4

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

molnam = "hirota"
@time qns, vals, vecs = westerenergies(parameters,0.5,nm,mc,0)
#testlines = (pred2lne(transitions))
#println(transitions)

#parameters .+= 0.5*(0.5 .- rand(size(parameters)))
#@time westerfit_handcoded()
