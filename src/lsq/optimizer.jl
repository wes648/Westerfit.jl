
function linereject(j,w,omc,unc,thres)
   filt = abs.(omc) .≤ (thres .* unc)
   jout = j[filt,:]
   wout = w[filt,filt]
   mout = omc[filt]
   return jout, wout, mout
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
   #@threads 
   for i in 1:length(perm)
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
   println("Initial RMS = $rms")
   goal = BLAS.nrm2(uncs)/√length(uncs)*ctrl["goal"]
   W = Diagonal(1.0 ./ uncs)
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
   #println(io,"")
   #println(io,"-------------------------------------")
   #println(io,"")
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
      #ρlm = lbmq_gain(βf,λlm,jtw*omc,omc,nomc)
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
         vecs .= nvecs
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
         iterationwriter(molnam,paramarray,rms,counter,λlm,βf,perm)
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
   covarmat = covarr(correl(H),puncs)
   #covarmat = inv(H)
   #params[1:15] .= paramrecov(params[1:15])
   #uncs[1:15] .= uncrecov(uncs[1:15],params[1:15])
   params[1:18], puncs[1:18] = fullrecov(params[1:18],puncs[1:18],ctrl["Irrep"])
   slλ = (@sprintf("%0.4f", log10(λlm)))
   outputfinal(molnam,ctrl,frms,counter,slλ,puncs,params,endpoint)
   if ctrl["overwrite"]==true
      println("Writing new input file at $molnam.inp. Previous file has moved to $molnam","1.inp")
      inpwriter(molnam, params, scales)
   end
   return params, covarmat, fomc, fcfrqs, vals
end
