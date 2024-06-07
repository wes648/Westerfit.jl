
function λgen(μ::Float64,er::Float64)::Float64
   ρ = 0.5
   er /= 1000.0 #convert err to GHz
   λ = ρ*μ*er #λF
   λ += (1.0 - ρ)*μ*er/(1+er) #λARC
   return λ
end

function jlister(inds)
   #finds all the unique J & σ pairs
   js = vcat(inds[:,1],inds[:,4])
   σs = vcat(inds[:,2],inds[:,5])
   temp = fill((0,0),size(js))
   for i in 1:size(js,1)
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

function rmscalc(vals,inds,ofreqs)
   cfreqs = zero(ofreqs)
   @threads for i in 1:size(cfreqs,1)
      cfreqs[i] = vals[inds[i,3],inds[i,2]+1] - vals[inds[i,6],inds[i,5]+1]
   end
   #println(cfreqs)
   omc = ofreqs - cfreqs
   rms = BLAS.nrm2(omc)/√length(omc)
   #rms = norm(omc)/length(omc)
   return rms, omc, cfreqs
end
#function quickfreqs(vals,inds,ofreqs)
#   cfreqs = zero(ofreqs)
#   @threads for i in 1:size(cfreqs,1)
#      cfreqs[i] = vals[inds[i,3],inds[i,2]+1] - vals[inds[i,6],inds[i,5]+1]
#   end
#   return cfreqs
#end
#construct jacobian
#function anaderiv(j,s,σ,vec,rp,rpid)
#   rp = zero(rp)
#   rp[rpid] = 1.0
#   U = ur(j,s,mcalc)
#   if σ==zero(σ)
#      U *= ut(mcalc,j,s)
#   end
#   mat = Matrix(U*Htsrmat(rp,j,s,mcalc,σ)*U)
#   out = transpose(vec)*mat*vec
#   return out
#end
function sumder(out,j,s,nf,rpid,prm,stg,ops,ms,qns)
   ind = rpid+1
   if ind ≤ length(stg)+15
      check = stg[ind-15]
      while check < zero(check)
         pm = prm[ind]
         #println(pm)
         #println(ops[:,ind-15])
         #out .+= tsrop(pm,ops[:,ind-15],j,s,nb,kb,mb,nk,kk,mk)
         out .+= tsr_op(1.0,j,s,qns,ms, ops[:, rpid-15] )
         ind += 1
         if ind-15 ≤ length(stg)
            check = stg[ind-15]
         else
            check = 0
         end
      end
   end
   return out
end

function derivmat_old(j,s,nf,rpid,prm,scl,stg,ops,nb,kb,mb,nk,kk,mk)
   #println(prm)
   if scl[rpid] < 0
   elseif rpid ≤ 4 #pure rot
      pr = zeros(4)
      pr[rpid] = 1.0
      out = hrsr(pr,zeros(4),zeros(3),j,s,nb,kb,nk,kk)
      out = kron(I(size(mk,1)),out)
   elseif 5 ≤ rpid ≤ 8 #spin-rot
      pr = zeros(4)
      pr[rpid-4] = 1.0
      out = hrsr(zeros(4),pr,zeros(3),j,s,nb,kb,nk,kk)
      out = kron(I(size(mk,1)),out)
   elseif 9 ≤ rpid ≤ 11 #qua
      pr = zeros(3)
      pr[rpid-8] = 1.0
      out = hrsr(zeros(4),zeros(4),pr,j,s,nb,kb,nk,kk)
      out = kron(I(size(mk,1)),out)
   elseif (rpid==12)||(rpid==14) # F or Vnf
      pr = zeros(3)
      pr[rpid-11] = 1.0
      out = kron(htorq(pr,nf,mb,mk), I(size(nk,1)))
   elseif rpid==13 # ρF
      out = kron(Diagonal(mk),Diagonal(kk))
   elseif rpid==15 # η 
      out = tsrop(1.0,0,0,0,0,1,1,0,0,j,s,nb,kb,mb,nk,kk,mk)
   else #user def
      out = tsrop(1.0,ops[:,rpid-15],j,s,nb,kb,mb,nk,kk,mk)
      out .= sumder(out,j,s,nf,rpid,prm,stg,ops,nb,kb,mb,nk,kk,mk)
   end
   return out
end
function derivmat(j,s,nf,rpid,prm,scl,stg,ops,ms,qns)
   if scl[rpid] < 0 #should this be ≤ 0 ???
   elseif rpid ≤ 4 #pure rot
      pr = zeros(4)
      pr[rpid] = 1.0
      out = hrot2(pr,qns)
      out = kron(I(length(ms)),out)
   elseif 5 ≤ rpid ≤ 9 #spin-rot
      pr = zeros(5)
      pr[rpid-4] = 1.0
      out = hsr(pr,j,s,qns)
      out = kron(I(length(ms)),out)
   elseif 10 ≤ rpid ≤ 12 #qua
      pr = zeros(3)
      pr[rpid-9] = 1.0
      out = hqu(pr,j,s,qns)
      out = kron(I(length(ms)),out)
   elseif rpid==13 # F
      pr = [1.;0.;0.;0.]
      out = kron(htor2(pr,ms), I(size(qns,1)))
   elseif rpid==16 # Vnf
      pr = [0.;0.;0.;1.]
      out = kron(htor2(pr,ms), I(size(qns,1)))
   elseif rpid==14 # ρzF
      out = kron(pa_op(ms,1), nz_op(qns,1))
   elseif rpid==15 # ρxF
      out = kron(pa_op(ms,1), npm_op(qns,1)) 
   elseif rpid==17 # ηz
      out = kron(pa_op(ms,1), sz_op(j,s,qns,1)) 
   elseif rpid==18 # ηx
      out = kron(pa_op(ms,1), spm_op(j,s,qns,1))
   else #user def
      #out = tsrop(1.0,ops[:,rpid-15],j,s,nb,kb,mb,nk,kk,mk)
      out = tsr_op(1.0,j,s,qns,ms,ops[:,rpid-18] )
      out .= sumder(out,j,s,nf,rpid,prm,stg,ops,ms,qns)
   end
   return out
end

function anaderiv(prm,scl,stg,rpid,ops,j,s,nf,ms,qns,vec)
   mat = derivmat(j,s,nf,rpid,prm,scl,stg,ops,ms,qns)
   out = transpose(vec)*mat*vec
   return diag(out)
end

function derivcalc(jlist,ops,ctrl,perm,vecs,prm,scl,stg)#removed nf call from here, fix in references
   #sd = Int(2*s+1)
   s = ctrl["S"]
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   #jmin = 0.5*iseven(Int(2*s+1))
   #println(jlist)
   #jmax = jlist[end,1]
   #jfd = Int(2*s+1)*Int(sum(2.0 .* collect(Float64,jmin:jmax) .+ 1.0))
   #vtd = Int(ctrl["vtmax"]+1)
   σcnt = σcount(nf)
   derivs = zeros(Float64,size(vecs,2),σcnt,length(perm))
   for sc in 1:σcnt
      #println(sc)
      σ = sc - 1
      msd = Int((2*mcalc+(σtype(nf,σ)==2)+1)*(2s+1))
      σt = σtype(nf,σ)
      #msd = Int(2*s+1)*mcd
      #mstrt, mstop = mslimit(nf,mcalc,σ)
      #jmsd = Int(msd*(2*jmax+1))
      #jsvd = Int(jfd*vtd)
      ms = msgen(nf,mcalc,σ)
      jsublist = jlist[isequal.(jlist[:,2],σ), 1] .* 0.5
      for j in jsublist
         #println(j)
         jd = Int(2.0*j) + 1
         sind, find = jvdest(j,s,ctrl["vtmax"]) 
         qns = qngen(j,s)
         vec = vecs[1:jd*msd,sind:find,sc]
         for i in 1:length(perm)
            pid = perm[i]
            ders = anaderiv(prm,scl,stg,pid,ops,j,s,nf,ms,qns,vec)
            derivs[sind:find,sc,i] = ders
         end#perm loop
      end #j loop
   end#σ loop
   return derivs
end#function

function derivcalc_all(ops,ctrl,perm,vecs,prm,scl,stg)
   s = ctrl["S"]
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   derivs = zeros(Float64,size(vecs,2),σcnt,length(perm))
   sd = Int(2.0*s+1.0)
   jmin = 0.5*iseven(sd)
   jmax = ctrl["Jmax"]
   jlist = collect(Float64,jmin:jmax)
   msd = Int((2*mcalc+(σtype(nf,σ)==2)+1)*(2s+1))
   for j in jlist
      jd = Int(2.0*j) + 1
      sind, find = jvdest(j,s,ctrl["vtmax"])
      qns = qngen(j,s)
      vec = veccs[1:jd*msd,sind:find,sc]
      for i in 1:length(perm)
         pid = perm[i]
         ders = anaderiv(prm,scl,stg,pid,ops,j,s,nf,ms,qns,vec)
         derivs[sind:find,sc,i] = ders
      end#perm loop
   end#j loop
   return derivs
end#function



function build_jcbn2!(jcbn,ops,jlist,inds,ctrl,vecs,params,perm,scals,stg)
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   #jlist = unique(vcat(inds[:,1:3],inds[:,4:6]))
   jcbn = zeros(Float64,size(inds,1),length(perm))
   deriv = derivcalc(jlist,ops,ctrl,perm,vecs,nf,params,scals,stg)
   @threads for p in 1:length(perm)
   @simd for a in 1:size(inds,1)
      jcbn[a,p] = deriv[inds[a,3],inds[a,2]+1,p] - deriv[inds[a,6],inds[a,5]+1,p]
   end
   end
   #@show jcbn
   return jcbn
end
function build_jcbn_sim(ops,inds,ctrl,vecs,params,perm,scals,stg)
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   jlist = #copy from tsrcalc
   #jlist = unique(vcat(inds[:,1:3],inds[:,4:6]))
   jcbn = zeros(Float64,size(inds)[1],length(perm))
   deriv = derivcalc(jlist,ops,ctrl,perm,vecs,nf,params,scals,stg)
   @threads for p in 1:length(perm)
   @simd for a in 1:size(inds,1)
      jcbn[a,p] = deriv[inds[a,3],inds[a,2]+1,p] - deriv[inds[a,6],inds[a,5]+1,p]
   end
   end
   return jcbn
end

function lbmq_gain(β,λ::Float64,g,omc,nomc)::Float64
   out = (0.5*transpose(β)*(λ*β + g))[1]
   out = (sum(omc .^2)-sum(nomc .^2)) / out
   return out
end

#function build_jcbn!(jcbn,ops,inds,s,ctrl,vecs,params,perm,scals)
#"""
#This builds the Jacobian based on the Hellmann–Feynman theorem.
#"""
#   nf = ctrl["NFOLD"]
#   mcalc = ctrl["mcalc"]
#   jcbn = zeros(Float64,size(inds)[1],length(perm))
#   @threads for a in 1:size(inds,1)
#      ju = 0.5*inds[a,1]
#      σu = inds[a,2]
#      nuk = ngen(ju,s)
#      kuk = kgen(ju,s)
#      muk = mgen(nf,mcalc,σu)
#      nub = permutedims(nuk)
#      kub = permutedims(kuk)
#      mub = permutedims(muk)
#      jl = 0.5*inds[a,4]
#      σl = inds[a,5]
#      nlk = ngen(jl,s)
#      klk = kgen(jl,s)
#      mlk = mgen(nf,mcalc,σl)
#      nlb = permutedims(nlk)
#      klb = permutedims(klk)
#      mlb = permutedims(mlk)
#      vecu = vecs[1:size(nuk,1)*size(muk,1),inds[a,3],σu+1]
#      vecl = vecs[1:size(nlk,1)*size(mlk,1),inds[a,6],σl+1]
#      @simd for i in 1:length(perm)
#         b = perm[i]
#         #dν/dOp = d/dOp (Eu - El)
#         jcbn[a,i]  = anaderiv(params,scals,b,ops,ju,s,nf,nub,kub,mub,nuk,kuk,muk,vecu)
#         jcbn[a,i] -= anaderiv(params,scals,b,ops,jl,s,nf,nlb,klb,mlb,nlk,klk,mlk,vecl)
#      end
#   end
#   return jcbn
#end
function tsrapprox(j,β)::Vector{Float64}
   for i in 1:length(β)
      j[:,i] .*= β[i]
   end
   #δ = sum(j,dims=2)
   return vec(sum(j,dims=2))
end

function linereject(j,w,omc,unc,thres)
   filt = abs.(omc) .≤ (thres .* unc)
   jout = j[filt,:]
   wout = w[filt,filt]
   mout = omc[filt]
   return jout, wout, mout
end

function build_hess!(hssn,dk,jtw,jcbn,weights)
   jtw = transpose(jcbn)*weights
   hssn = jtw*jcbn
   #Threads.@threads for i in 1:size(hssn)[1]
   #   dk[i,i] = norm(hssn[:,i])
   #end
   return hssn, dk, jtw
end
function build_hess(jtw,jcbn,weights)
   jtw = transpose(jcbn)*weights
   hssn = jtw*jcbn
   #Threads.@threads for i in 1:size(hssn)[1]
   #   dk[i,i] = norm(hssn[:,i])
   #end
   #println(hssn)
   return hssn, jtw
end

function approx2dirdrv!(K,β,jcbn,weights,nlist,inds,params,perm,omc,ofreqs)
   h = 0.1
   params[perm] += h*β
   #pvals,pvecs = limeigcalc(nlist, inds, params)
   pvals,pvecs, = tsrcalc(ctrl,prm,stg,cdo,nf,vtm,mcalc,jlist,s,sd,σ)
   prms, pomc, = rmscalc(pvals, inds, ofreqs)
   K = (2/h)*((pomc-omc)/h - jcbn*β)
   return K
end
function lbmq_acc!(K,β,jcbn,weights,nlist,inds,params,perm,omc,ofreqs,λ)
   jtw = transpose(jcbn)*weights
   jtj = jtw*jcbn
   A = Hermitian(jtj + λ*Diagonal(jtj))#transpose(dk)*dk
   while isposdef(A)==false #this could be tidier
      λ = max(2.0*λ,1.0E-24)
      A .= Hermitian(jtj + λ*Diagonal(jtj))
   end
   A = cholesky!(Hermitian(A))
   K = approx2dirdrv!(K,β,jcbn,weights,nlist,inds,params,perm,omc,ofreqs)
   X = jtw*K
   β2 = zero(β)
   β2 = ldiv!(β2, A, -X)
   if (2*norm(β2)/norm(β))>0.5
      β2 *= 0.5*0.5*norm(β)/norm(β2)
      #β2 .*= 0.0
   end
   β += 0.5*β2
   return β
end

function lbmq_step!(β,H,grad, λ)
#   jtw = transpose(jcbn)*weights
#   jtj = jtw*jcbn
   A = Hermitian(H + λ*Diagonal(H))#transpose(dk)*dk
   while isposdef(A)==false #this could be tidier
      λ = max(2.0*λ,1.0E-24)
      A .= Hermitian(H + λ*Diagonal(H))
   end
   A = cholesky!(A)
   β .= ldiv!(β, A, -grad)
   return β,λ
end

function wellcon_model(t,p)
   @. p[1]*exp(p[2]*t)
end

function wellcon_acc()
#well conditioned accelerator. If λ is very small and previous 4 iterations
# have same sign, does an exponential fit and extrapolates to the end of the fit
   xs = collect(1:count)
   p0 = zeros(2)
   γ = zeros(length(perm))
   @threads for i in 1:length(perm)
      y = βset[i,:]
      if allequal(sign.(y))
         p0[1] = y[1]
         fit = curve_fit(wellcon_model,xs,y,p0)
         a = fit.param[1]
         k = fit.param[2]
         γ[i] = -a/(exp(k)-1.0) - sum(y)
      else
         γ[i] = 0.0
      end
   end
   return γ
end

function turducken_acc(λ::Float64,β::Array{Float64},h)::Float64
   βt = transpose(β)
   out = 1.0 + λ*βt*β/(βt*h*β)
   return out
end

function harshfilt!(β,param,scals)
   tol = 5.0e-1
   for i in 1:length(param)
      test = param[i]*scals[i]*tol
      if abs(β[i]) > test
         β[i] = sign(β[i])*test
      end
   end
   return β
end
function trfilter!(β,h,Δ)
   nrm = BLAS.nrm2(β .* diag(h))
   if nrm > Δ
      β .*= Δ/nrm
   end
   return β
end

#function lbmq_turducken!(βf,D,H,jtw,omc,λ,nlist,inds,nparams,perm,ofreqs)
function lbmq_turducken!(H,J,jtw,omc,λ,Δ,nlist,inds,nparams,scls,perm,ofreqs,rms,stg,cdo,ctrl)
   tdncount = ctrl["turducken"]
   A = Hermitian(H + λ*Diagonal(H))
   while isposdef(A)==false #this could be tidier
      λ = max(2.0*λ,1.0E-24)
      #println(λ)
      A = Hermitian(H + λ*Diagonal(H))
   end
   if isinf(λ)
      @warn "LB-MQ Matrix not pos-def!"
      println("Make sure you aren't trying to optimize a parameter with value of 0.0.")
      println("This code is about to crash")
   end
   A = cholesky!(A)
   β = zeros(Float64,length(perm),tdncount)
   β[:,1] .= ldiv!(β[:,1], A, jtw*omc) .* scls[perm]
   for i in 2:tdncount
      nparams[perm] .+= β[:,i-1]
      vals,nvecs, = tsrcalc2(nparams,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
      nrms, omc, = rmscalc(vals,inds,ofreqs)
      β[:,i] .= ldiv!(β[:,i], A, jtw*omc) .* scls[perm]
   end
   βf = sum(β,dims=2)
   nparams[perm] .+= β[:,end]
   vals, nvecs = tsrcalc2(nparams,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
   nrms, omc, = rmscalc(vals,inds,ofreqs)
   return βf,λ,omc,nrms,vals,nvecs, nparams
end

function aitkenδ(γ)#this was a weird accelerator idea
   @. (γ[:,end-2]*γ[:,end] - γ[:,end-1]^2)/(γ[:,end-2]+γ[:,end] - 2.0*γ[:,end-1])
end
function paramunc(H,W,perm,omc)
   uncs = diag(inv(Symmetric(H)))
   uncs .*= (omc' * W * omc)/(length(omc)-length(perm))
   return □rt.(uncs)
end
function permdeterm(scls,stgs)
   out = collect(1:length(scls))[(scls .> 0) .* (vcat(ones(18),stgs) .> 0)]
end

function lbmq_opttr(ctrl,nlist,ofreqs,uncs,inds,params,scales,cdo,stg,molnam)
   #vals,vecs = limeigcalc(nlist, inds, params)
   #S = ctrl["S"]
   #println(inds)
   #sd = Int(2*S+1)
   vals,vecs, = tsrcalc2(params,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
   LIMIT = ctrl["maxiter"]


   paramarray = zeros(Float64, length(params), LIMIT+1)
   paramarray[:,1]=params
   oparams = paramarray[:,1]

   rms, omc, = rmscalc(vals, inds, ofreqs)
   #perm,n = findnz(sparse(scales))
   perm = permdeterm(scales,stg)
   println(perm)
   #println(params)
   #println(omc)
   #println(nlist)
   println("Initial RMS = $rms")
   goal = sum(uncs)/length(uncs)*0.00000
   W = diagm(0=>(uncs .^ -1))
   #RHOTHRES = -1.0E-6
   ϵ0 = 0.1E-6
   ϵ1 = 0.1E-6
   μlm = ctrl["λlm0"]#(rms + rms^2)#*0.0
   λlm = λgen(μlm, rms) 
   oλlm = λlm
   println("Initial λ = $λlm")
   Δlm = nrm2(params[perm])/length(perm)
   counter = 0
   BAD = 0

   io = open("$molnam.out", "a")
   println(io,"Initial RMS = $rms MHz")
   println(io,"Initial λ = $λlm")
   println(io,"")
   println(io,"-------------------------------------")
   println(io,"")
   close(io)

   nparams = copy(params)
   #puncs = zero(perm)
   βf = zero(perm) #step
   J = zeros(Float64,size(inds,1),length(perm)) #Jacobian
   jtw = zero(omc) #jtwient
   J = build_jcbn2!(J,cdo,nlist,inds,ctrl,vecs,params,perm,scales,stg)
   #J, w, omc = linereject(J,W,omc,uncs,ctrl["REJECT"])
   H, jtw = build_hess(jtw,J,W)
   #println(H)
   if true ∈ isnan.(H)
      println("FUCKING FUCKING FUCK. NaN in Hessian")
   end
   #uncs = paramunc(uncs,H,perm,omc)
   endpoint = "not yet"
   converged=false
   while converged==false
      #=
      β,λlm = lbmq_step!(β,H,grad,λlm)
      if true
         β .= lbmq_acc!(K,β,J,W,nlist,inds,params,perm,omc,ofreqs,λlm)
      end
      #if (norm((H^(-1/2))*β)>Δlm)&&(λ!=0.0)
      #   β *= Δlm/norm(D*β)
      #end
      nparams[perm] .= params[perm] + β #+ (inv(H)^2)*β
      vals, nvecs = limeigcalc(nlist, inds, nparams)
      nrms, nomc = rmscalc(vals,inds,ofreqs)=#
      #oparams = copy(params)
      λlm = λgen(μlm, rms) 
      βf,λlm,nomc,nrms,vals,nvecs,nparams = lbmq_turducken!(H,J,
         jtw,omc,λlm,Δlm,nlist,inds,copy(params),scales,perm,ofreqs,rms,stg,cdo,ctrl)
      check = abs(nrms-rms)/rms
      #println(βf)
      ρlm = lbmq_gain(βf,λlm,jtw*omc,omc,nomc)
   #println()
   #println(ρlm)
   #println()
      if nrms < rms#*(0.95 + 0.3*exp(-0.6*BAD))
         if nrms < rms
            BAD = max(0,BAD-1)
         else
            println("This might be a bad step")
            BAD += 1
         end
   #    if ρlm > 0.0#-1.0e-7 #
         #println(βf)
	 #μlm *= (nrms/rms)^2
         rms = nrms
         omc = nomc
         params .= nparams
         #println(params)
         #vecs .= nvecs
         #@time J = build_jcbn!(J,cdo,inds,S,ctrl,vecs,params,perm,scales)
         @time J = build_jcbn2!(J,cdo,nlist,inds,ctrl,vecs,params,perm,scales,stg)
         #J, w = linereject(J,W,omc,uncs,ctrl["REJECT"])
         H, jtw = build_hess(jtw,J,W)
         #println(diag(H))
         counter += 1
         paramarray[:,counter+1] = params
         #sρlm = (@sprintf("%0.4f", ρlm))
         srms = (@sprintf("%0.4f", rms))
         slλ = (@sprintf("%0.4f", log10(λlm)))
         #sΔ = (@sprintf("%0.6f", Δlm))
         scounter = lpad(counter,3)
         println("After $scounter iterations, RMS = $srms, log₁₀(λ) = $slλ")
         iterationwriter(molnam,paramarray,srms,scounter,slλ,βf,perm)
         #println(H^(-1/2))
         #println(params)
         #println(diag(H))
         #λlm = λgen(μlm, rms) 
         #=if λlm < oλlm
            μlm /= 30.0
         else
            μlm /= 3.0 #20.0
         end=#
         μlm /= 30.0
         #oλlm = λlm
         #Δlm *= 6.0
      else
         #params .= oparams
         μlm = max(4.0*μlm,1.0E-24)
         #Δlm = max(0.90*Δlm,0.0001)
      end
      #ρlm = lbmq_gain(β,λlm,jtw*omc,rms,nrms)
      #if ρlm ≥ 0.75
      #   Δlm *= 2.0
      #else
      #   #Δlm = max(0.90*Δlm,0.0001)
      #end
      if (rms ≤ goal)#&&(counter > 1)
         println("A miracle has come to pass. The fit has converged")
         endpoint = "converge"
         break
      elseif (check < ϵ0)
         println("The RMS has stopped decreasing. Hopefully it is low")
         #uncs = paramunc(uncs,H,perm,omc)
         #println(omc)
         #println(uncs)
         endpoint = "RMS"
         break
      elseif (norm(βf))<ϵ1*(norm(params[perm])+ϵ1)
         slλ = (@sprintf("%0.4f", log10(λlm)))
         println("It would appear step size has converged. log₁₀(λ) = $slλ")
         #uncs = paramunc(uncs,H,perm,omc)
         #println(uncs)
         endpoint = "step size"
         break
      elseif (λlm > 1.0e+9)&&(Δlm == 0.0)
         println("λlm exceeded threshold.")
         println("If you were using the turducken, try again without it")
         endpoint = "LMthresh"
         break
      elseif counter ≥ LIMIT
         println("Alas, the iteration count has exceeded the limit")
         #println(omc)
         endpoint = "iter"
         break
      else
      end #check if
   end#while
   frms, fomc, fcfrqs = rmscalc(vals, inds, ofreqs)
   puncs = zeros(size(params))
   puncs[perm] = paramunc(H,W,perm,omc)
   puncs_forsim = copy(puncs)
   #params[1:15] .= paramrecov(params[1:15])
   #uncs[1:15] .= uncrecov(uncs[1:15],params[1:15])
   params[1:18], puncs[1:18] = fullrecov(params[1:18],puncs[1:18],ctrl["Irrep"])
   slλ = (@sprintf("%0.4f", log10(λlm)))
   outputfinal(molnam,ctrl,frms,counter,slλ,puncs,params,endpoint)
   if ctrl["overwrite"]==true
      println("Writing new input file at $molnam.inp. Previous file has moved to $molnam","1.inp")
      inpwriter(molnam, params)
   end
   return params, puncs_forsim, fomc, fcfrqs, vals
end
