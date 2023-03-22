
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

function rmscalc(vals,inds,ofreqs)
   cfreqs = zero(ofreqs)
   @threads for i in 1:size(cfreqs)[1]
      cfreqs[i] = vals[inds[i,3],inds[i,2]+1] - vals[inds[i,6],inds[i,5]+1]
   end
   #println(cfreqs)
   omc = ofreqs - cfreqs
   rms = BLAS.nrm2(omc)/√length(omc)
   #rms = norm(omc)/length(omc)
   return rms, omc
end
#construct jacobian
function anaderiv(j,s,σ,vec,rp,rpid)
   rp = zero(rp)
   rp[rpid] = 1.0
   U = ur(j,s,mcalc)
   if σ==zero(σ)
      U *= ut(mcalc,j,s)
   end
   mat = Matrix(U*Htsrmat(rp,j,s,mcalc,σ)*U)
   out = transpose(vec)*mat*vec
   return out
end
function sumder(out,j,s,nf,rpid,prm,scl,ops,nb,kb,mb,nk,kk,mk)
   ind = rpid+1
   if ind ≤ length(scl)
      check = scl[ind]
      while check < zero(check)
         pm = prm[ind]
         out .+= tsrop(pm,ops[:,ind-15],j,s,nb,kb,mb,nk,kk,mk)
         ind += 1
         check = scl[ind]
      end
   end
   return out
end

function derivmat(j,s,nf,rpid,prm,scl,ops,nb,kb,mb,nk,kk,mk)
   if rpid ≤ 4 #pure rot
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
      out .= sumder(out,j,s,nf,rpid,prm,scl,ops,nb,kb,mb,nk,kk,mk)
   end
   return out
end
function anaderiv(prm,scl,rpid,ops,j,s,nf,nb,kb,mb,nk,kk,mk,vec)
   mat = derivmat(j,s,nf,rpid,prm,scl,ops,nb,kb,mb,nk,kk,mk)
   out = transpose(vec)*mat*vec
   return out
end

function lbmq_gain(β,λ,g,rms,nrms)
   out = 0.5*transpose(β)*(λ*β-g)
   out = (rms-nrms)/out
   return out
end

function build_jcbn!(jcbn,ops,inds,s,ctrl,vecs,params,perm,scals)
"""
This builds the Jacobian based on the Hellmann–Feynman theorem.
"""
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   jcbn = zero(jcbn)
   @threads for a in 1:size(inds,1)
      ju = 0.5*inds[a,1]
      σu = inds[a,2]
      nuk = ngen(ju,s)
      kuk = kgen(ju,s)
      muk = mgen(nf,mcalc,σu)
      nub = Matrix(transpose(nuk))
      kub = Matrix(transpose(kuk))
      mub = Matrix(transpose(muk))
      jl = 0.5*inds[a,4]
      σl = inds[a,5]
      nlk = ngen(jl,s)
      klk = kgen(jl,s)
      mlk = mgen(nf,mcalc,σl)
      nlb = Matrix(transpose(nlk))
      klb = Matrix(transpose(klk))
      mlb = Matrix(transpose(mlk))
      vecu = vecs[1:size(nuk,1)*size(muk,1),inds[a,3],σu+1]
      vecl = vecs[1:size(nlk,1)*size(mlk,1),inds[a,6],σl+1]
      @simd for i in 1:length(perm)
         b = perm[i]
         #dν/dOp = d/dOp (Eu - El)
         jcbn[a,i]  = anaderiv(params,scals,b,ops,ju,s,nf,nub,kub,mub,nuk,kuk,muk,vecu)
         jcbn[a,i] -= anaderiv(params,scals,b,ops,jl,s,nf,nlb,klb,mlb,nlk,klk,mlk,vecl)
      end
   end
   return jcbn
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
   pvals,pvecs, = tsrcalc(prm,stg,cdo,nf,vtm,mcalc,jlist,s,sd,σ)
   prms, pomc = rmscalc(pvals, inds, ofreqs)
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

#function lbmq_turducken!(βf,D,H,jtw,omc,λ,nlist,inds,nparams,perm,ofreqs)
function lbmq_turducken!(H,jtw,omc,λ,nlist,inds,nparams,perm,ofreqs,rms,stg,cdo,ctrl)
   if rms > 3000.
      tdncount = 1
   else
      tdncount = 3
   end
   A = Hermitian(H + λ*Diagonal(H))
   while isposdef(A)==false #this could be tidier
      λ = max(2.0*λ,1.0E-24)
      println(λ)
      A = Hermitian(H + λ*Diagonal(H))
   end
   if isinf(λ)
      println("LVMB Matrix not pos-def!") 
      println("Make sure you aren't trying to optimize a parameter with value of 0.0.")
      println("This code is about to crash")
   end
   A = cholesky!(A)
   β = zeros(Float64,length(perm),tdncount)
   β[:,1] .= ldiv!(β[:,1], A, jtw*omc)
   for i in 2:tdncount
      nparams[perm] .+= β[:,i-1]
      #vals, nvecs = limeigcalc(nlist, inds, nparams)
      vals,nvecs, = tsrcalc2(nparams,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
      nrms, omc = rmscalc(vals,inds,ofreqs)
      β[:,i] .= ldiv!(β[:,i], A, jtw*omc)
   end
   #δβ = aitkenδ(β)
   βf = sum(β,dims=2) #.+ δβ
   nparams[perm] .+= β[:,end] #.+ δβ
   vals, nvecs = tsrcalc2(nparams,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
   nrms, omc = rmscalc(vals,inds,ofreqs)
   return βf,λ,omc,nrms,vals,nvecs, nparams
end

function aitkenδ(γ)#the fuck is this like literally what is this for
   @. (γ[:,end-2]*γ[:,end] - γ[:,end-1]^2)/(γ[:,end-2]+γ[:,end] - 2.0*γ[:,end-1])
end

function paramunc!(uncs,H,perm,omc)
   uncs = □rt.(diag(inv(Symmetric(H)))) .* (sum(omc .^2))/(length(omc)-length(perm))
   return uncs
end

function lbmq_opttr(ctrl,nlist,ofreqs,uncs,inds,params,scales,cdo,stg)
   #vals,vecs = limeigcalc(nlist, inds, params)
   S = ctrl["S"]
   #println(inds)
   sd = Int(2*S+1)
   vals,vecs, = tsrcalc2(params,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
   oparams = params
   rms, omc = rmscalc(vals, inds, ofreqs)
   perm,n = findnz(sparse(scales))
   println(perm)
   println(params)
   println(omc)
   println("Initial RMS = $rms")
   goal = sum(uncs)/length(uncs)*0.00000
   W = diagm(0=>(uncs .^ -1))
   #RHOTHRES = -1.0E-6
   ϵ0 = 0.1E-12
   ϵ1 = 0.1E-16
   LIMIT = 50
   μlm = rms + rms^2
   λlm = μlm*rms/(1.0 + rms)
   Δlm = 1.0E+2
   Δlm *= length(perm)
   counter = 0
   stoit = 5
   #println(omc)
   nparams = copy(params)
   uncs = zero(perm)
   βf = zero(perm) #step
   J = zeros(Float64,size(inds)[1],length(perm)) #Jacobian
   jtw = zero(omc) #jtwient
   K = zero(omc) #acceleration correction
   #H = zeros(Float64,length(perm),length(perm)) #Hessian
   #D = zero(H) #Elliptical Trust Region
   J = build_jcbn!(J,cdo,inds,S,ctrl,vecs,params,perm,scales)
   H, jtw = build_hess(jtw,J,W)
   if true ∈ isnan.(H)
      println("FUCKING FUCKING FUCK")
   end
   #uncs = paramunc!(uncs,H,perm,omc)
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
      λlm = μlm*rms/(1.0 + rms)
   βf,λlm,nomc,nrms,vals,nvecs,nparams = lbmq_turducken!(H,
                  jtw,omc,λlm,nlist,inds,copy(params),perm,ofreqs,rms,stg,cdo,ctrl)
      check = abs(nrms-rms)/rms
      if nrms < rms
         println(βf)
         rms = nrms
         omc .= nomc
         params .= nparams
         println(params)
         #vecs .= nvecs
         J = build_jcbn!(J,cdo,inds,S,ctrl,vecs,params,perm,scales)
         H, jtw = build_hess(jtw,J,W)
         counter += 1
         srms = (@sprintf("%0.4f", rms))
         slλ = (@sprintf("%0.4f", log10(λlm)))
         #sΔ = (@sprintf("%0.6f", Δlm))
         scounter = lpad(counter,3)
         println("After $scounter interations, RMS = $srms, log₁₀(λ) = $slλ")#, Δₖ = $sΔ")
         #println(H^(-1/2))
         #println(params[perm])
         #λlm = 0.0
         μlm /= 20.0
         stoit = 0
      else
         #params .= oparams
         μlm = max(4.0*μlm,1.0E-24)
      end
      #ρlm = lbmq_gain(β,λlm,jtw*omc,rms,nrms)
      #if ρlm ≥ 0.75
      #   Δlm *= 2.0
      #else
      #   #Δlm = max(0.90*Δlm,0.0001)
      #end
      if (rms ≤ goal)#&&(counter > 1)
         println("A miracle has come to pass. The fit has converged")
         break
      elseif (check < ϵ0)
         println("The RMS has stopped decreasing. Hopefully it is low")
         uncs = paramunc!(uncs,H,perm,omc)
         println(uncs)
         break
      elseif (norm(βf))<ϵ1*(norm(params[perm])+ϵ1)
         slλ = (@sprintf("%0.4f", log10(λlm)))
         println("It would appear step size has converged. log₁₀(λ) = $slλ")
         uncs = paramunc!(uncs,H,perm,omc)
         println(uncs)
         break
      elseif counter ≥ LIMIT
         println("Alas, the iteration count has exceeded the limit")
         #println(omc)
         break
      else
         #write update to file
      end #check if
   end#while
   return params, vals
end
