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
using LinearAlgebra
using LinearAlgebra.BLAS
using LinearAlgebra.LAPACK
using LineSearches
using SparseArrays
include("/home/wes/files/westerfit/src/common.jl")
include("/home/wes/files/westerfit/src/hamiltonian.jl")
include("/home/wes/files/westerfit/src/intensities.jl")
include("/home/wes/files/westerfit/src/WIGXJPF.jl")
using .WIGXJPF

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


function assign(nmax,j,s,σ,mcalc,mmax,vals,vecs)
   #determine which indices we want to assign
   nlist = Δlist(j,s)
   ndegns = @. 2*nlist + 1
   mcind = mcalc+1
   mcd = 2*mcalc+1
   offset = zeros(Int64,length(nlist))
   offset[2:end] = (mcd) .* ndegns[1:end-1]
   goals = zeros(Int,0,1)
   off = 0
   for n in nlist
      start = off + Int((mcalc - mmax)*(2*n+1)) + 1
      finsh = off + Int((mcalc + mmax + 1)*(2*n+1))
      part = collect(start:finsh)
      goals = vcat(goals,part)
      off += Int((2*mcalc+1)*(2*n+1))
   end
   assigned = leadfact(vals,copy(vecs))
   vals = vals[assigned]
   vecs = vecs[:,assigned]
   vals = vals[goals]
   vecs = vecs[:,goals]
   QNs = qngen(j,s,mmax,σ)
   return QNs, vals, vecs
end

function leadfact(vals,vecs)
   vecs = abs.(vecs)
   c = zeros(Int,length(vals))#,2)
   for i in 1:length(vals)
      t = argmax(vecs)
      s = t[1] #state number that this eigenvalue maps to
      k = t[2] #eigenvalue number from energetic sorting
      vecs[s,:] = zeros(Float64,size(vecs[1,:])) #prevents state from being double assigned
      vecs[:,k] = zeros(Float64,size(vecs[:,2]))
      c[k] = s
#      c[k,2] = i
   end
   perm = sortperm(c)
   return perm
end

function tsrdiag(pr,j,s,nmax,mcalc,mmax,σ)
   #U = ur(j,s,mcalc)
   if σ==0||σ==0.0
      U = ur(j,s,mcalc)*ut(mcalc,j,s)
      H = Matrix(U*Htsr(pr,j,s,mcalc,σ)*U)
   else
      H = Matrix(Htsr(pr,j,s,mcalc,σ))
   end
   vals, vecs = LAPACK.syev!('V', 'U', H)
   qns, vals, vecs = assign(nmax,j,s,σ,mcalc,mmax,vals,vecs)
   return qns, vals, vecs
end

function tsrcalc(prm,s,nmax,mcalc,mmax,σ)
   #prm = PAM2RAM(prm)
   prm = paramprep(prm)
   #prm = paramshift(inprms)
   if isodd(Int64(2*s+1))
      jmin = 0.0
   else
      jmin = 0.5
   end
   jmax = nmax - s
   jarray = collect(Float64,jmin:jmax)
   jd = Int((2*s+1)*sum(2 .* jarray .+ 1))
   outvals = zeros(Float64,Int(jd*(2*mmax+1)),2*mmax+1)
   outvecs = zeros(Float64,Int((2*s+1)*(2*jmax+2)*(2*mcalc+1)),Int(jd*(2*mmax+1)),2*mmax+1)
   outqns = zeros(Float64,Int(jd*(2*mmax+1)),6)
   for i in 1:length(jarray)
      j = jarray[i]
      sind, find = jinds(j,s)
      tqns, tvals, tvecs = tsrdiag(prm,j,s,nmax,mcalc,mmax,σ)
      outvals[sind:find] = tvals
      outvecs[1:Int((2*s+1)*(2*j+1)*(2*mcalc+1)),sind:find] = tvecs
      outqns[sind:find,:] = tqns
   end
   return outqns, outvals, outvecs
end

function westersim(prm,s,nmax,mcalc,mmax)
   @time amq, ama, amv = tsrcalc(parameters,0.5,nm,mc,0,0)
   @time atransitions = tracalc(nm,0.5,mc,0,amq,ama,amv)
   println("E states!")
   @time emq, ema, emv = tsrcalc(parameters,0.5,nm,mc,0,1)
   qns  = cat(amq,emq, dims=3)
   vals = cat(ama,ema, dims=3)
   vecs = cat(amv,emv, dims=3)
   @time etransitions = tracalc(nm,0.5,mc,1,emq,ema,emv)
   transitions = vcat(atransitions,etransitions)
   return transitions, qns, vals, vecs
end

################################################################################
############                  westerfit                  ############
################################################################################

function lineprep(lns)
   #converts the input file into a more code friendly format
   #           1  2  3   4   5  6  7  8  9   10 11 12  13   14
   #input  = [ju nu kau kcu mu σu jl nl kal kcl ml σl freq unc]
   #           1  2   3   4  5   6
   #output = [ju σu indu jl σl indl]
   qunus = lns[:,1:12]
   freqs = lns[:,13]
   uncs = lns[:,end]
   inds = zeros(Int64,size(lns)[1],6)
   inds[:,1] = Int64.(2 .* qunus[:,1])
   inds[:,2] = Int64.(qunus[:,6])
   inds[:,3] = qn2ind.(qunus[:,1],0.5,qunus[:,2],qunus[:,3],qunus[:,4])
   inds[:,4] = Int64.(2 .* qunus[:,7])
   inds[:,5] = Int64.(qunus[:,12])
   inds[:,6] = qn2ind.(qunus[:,7],0.5,qunus[:,8],qunus[:,9],qunus[:,10])
   #inds = vcat(inds[:,1:2], inds[:,3:4])
   return inds, freqs, uncs
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
   rp = paramprep(rp)
   jmax = 0.5*maximum(jlist[:,1])
   #println(jlist)
   σs = σcount(NFOLD)
   jd = (2*S+1)*sum(2 .* collect((0.5*isodd(2*S)):jmax) .+ 1)
   evals = zeros(Float64,Int(jd*(2*mmax+1)),σs)
   evecs = zeros(Float64,Int((2*S+1)*(2*jmax+1)*(2*mcalc+1)),Int(jd*(2*mmax+1)),σs)
   eqns = zeros(Float64,Int(jd*(2*mmax+1)),6,σs)
   for i in 1:size(jlist)[1]
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
   for a in 1:size(inds)[1]
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
   U = ur(j,s,mcalc)
   if σ==0
      U *= ut(mcalc,j,s)
   end
   mat = Matrix(U*Htsr(drp,j,s,mcalc,σ)*U)
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
   U = ur(j,s,mcalc)
   if σ==0
      U *= ut(mcalc,j,s)
   end
   mat = Matrix(U*Htsr(drp,j,s,mcalc,σ)*U)
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
   A = Hermitian(jtw*jcbn + lbmq_prm*I)
   A = cholesky(A)
   X = jtw*omc
   β = ldiv!(β, A, X)
   return β
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
   LIMIT = 200
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
   λ = 1.0
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



c = 29979.2458 #MHz to cm-1 converter
F = 292574.52
V3 = 4195821.0
A = 82756.97
B = 10145.212
C = 9426.952
δ = 0.45
#Dab = 20.0
ϵzz = 254.3618
ϵxx = 385.4886
ϵyy = -415.9114
ϵxz = 0.0
ΔN = 0.01188
η = 0.0
μa = 1.0
μb = 0.5

#parameters = [ A;   B;   C;   δ;   F;   V3; ϵzz; ϵxx; ϵyy; ϵxz;  η]
#trsscales =  [0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01]
parameters = [83316.10066771678, 14497.783454870654, 4733.339691158756,
      -1.1414038835204492, 292650.6595414094, 4.195803378991427e6,
      306.43192861068763, 297.11625650644,
      -22.797297784556505, 0.0, 0.0, 33.90230680364543]

trsscales = fill(1.0,12)

"""
testlines = [1.5  1  0  1  0  0  0.5  0  0  0  0  0  19419.0217  0.4;
 0.5  1  0  1  0  0  0.5  0  0  0  0  0  19450.6917  0.4;
 2.5  2  0  2  0  0  0.5  1  0  1  0  0  38819.1124  0.4;
 2.5  2  0  2  0  0  1.5  1  0  1  0  0  38847.8178  0.4;
 1.5  2  0  2  0  0  0.5  1  0  1  0  0  38865.0721  0.4;
 1.5  2  0  2  0  0  1.5  1  0  1  0  0  38894.2759  0.4;
 2.5  3  0  3  0  0  1.5  2  0  2  0  0  58287.5536  0.4;
 3.5  3  0  3  0  0  2.5  2  0  2  0  0  58272.0583  0.4;
 1.5  1  0  1  1  1  0.5  0  0  0  1  1  19104.1433  0.4;
 0.5  1  0  1  1  1  0.5  0  0  0  1  1  19137.5883  0.4;
 2.5  2  0  2  1  1  1.5  1  0  1  1  1  38250.8564  0.4;
 1.5  2  0  2  1  1  0.5  1  0  1  1  1  38263.2012  0.4;
 1.5  2  0  2  1  1  1.5  1  0  1  1  1  38295.7924  0.4;
 2.5  2  0  2  1  1  0.5  1  0  1  1  1  38228.7521  0.4]
"""

INTTHRESHOLD = .000001
#=
parameters:
     1  2  3   4   5  6   7   8  9 10 11 12
prs =[A; B; C; Dab; Fr; ρ; V3; ao; a; b; d; η]

λs.vectors[abs.(λs.vectors) .< 1E-9] .= 0.0
=#
mc = 10
nm = 3

@time transitions, qns, vals, vecs = westersim(parameters,0.5,nm,mc,0)
testlines = abs.(pred2lne(transitions))

parameters .+= 0.5*(0.5 .- rand(size(parameters)))
@time westerfit_optim()
