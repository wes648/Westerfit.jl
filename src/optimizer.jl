
function rmscalc(vals,inds,ofreqs)
   cfreqs = zero(ofreqs)
   for i in 1:size(cfreqs)[1]
      cfreqs[i] = vals[inds[i,3],inds[i,2]+1] - vals[inds[i,6],inds[i,5]+1]
   end
   #println(cfreqs)
   omc = ofreqs - cfreqs
   rms = √(sum(omc .^ 2)/length(omc))
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
   mat = Matrix(U*Htsr(rp,j,s,mcalc,σ)*U)
   out = transpose(vec)*mat*vec
   return out
end
function lbmq_gain(β,λ,g,rms,nrms)
   out = 0.5*transpose(β)*(λ*β-g)
   out = (rms-nrms)/out
   return out
end

function build_jcbn!(jcbn,inds,vecs,params,perm)
"""
This builds the Jacobian based on the Hellmann–Feynman theorem.
"""
   @threads for a in 1:size(inds)[1]
      ju = 0.5*inds[a,1]
      jl = 0.5*inds[a,4]
      σu = inds[a,2]
      σl = inds[a,5]
      vecu = vecs[1:Int((2*S+1)*(2*ju+1)*(2*mcalc+1)),inds[a,3],σu+1]
      vecl = vecs[1:Int((2*S+1)*(2*jl+1)*(2*mcalc+1)),inds[a,6],σl+1]
      for i in 1:length(perm)
         b = perm[i]
         jcbn[a,i] = anaderiv(ju,S,σu,vecu,params,b) - anaderiv(jl,S,σl,vecl,params,b)
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
   return hssn, jtw
end

function approx2dirdrv!(K,β,jcbn,weights,nlist,inds,params,perm,omc,ofreqs)
   h = 0.1
   params[perm] += h*β
   pvals,pvecs = limeigcalc(nlist, inds, params)
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

#function lbmq_turducken!(βf,D,H,jtw,omc,λ,nlist,inds,nparams,perm,ofreqs)
function lbmq_turducken!(H,jtw,omc,λ,nlist,inds,nparams,perm,ofreqs)
   A = Hermitian(H + λ*Diagonal(H))
   while isposdef(A)==false #this could be tidier
      λ = max(2.0*λ,1.0E-24)
      A .= Hermitian(H + λ*Diagonal(H))
   end
   A = cholesky!(A)
   β = zeros(Float64,size(perm))
   β .= ldiv!(β, A, jtw*omc)
   βf = copy(β)
   nparams[perm] .+= β
   vals, nvecs = limeigcalc(nlist, inds, nparams)
   nrms, nomc = rmscalc(vals,inds,ofreqs)
   β .= ldiv!(β, A, jtw*nomc)
   βf .+= β
   nparams[perm] .+= β
   vals, nvecs = limeigcalc(nlist, inds, nparams)
   nrms, nomc = rmscalc(vals,inds,ofreqs)
   β .= ldiv!(β, A, jtw*nomc)
   βf .+= β
   nparams[perm] .+= β
   vals, nvecs = limeigcalc(nlist, inds, nparams)
   nrms, nomc = rmscalc(vals,inds,ofreqs)
   return βf,λ,nomc,nrms,vals,nvecs, nparams
end

function paramunc!(uncs,H)
   uncs = diag(inv(H))
   return uncs
end

function lbmq_opttr(nlist,ofreqs,uncs,inds,params,scales,λ)
   vals,vecs = limeigcalc(nlist, inds, params)
   oparams = params
   rms, omc = rmscalc(vals, inds, ofreqs)
   perm,n = findnz(sparse(scales))
   println(perm)
   println("Initial RMS = $rms")
   goal = sum(uncs)/length(uncs)*0.00000
   W = diagm(0=>(uncs .^ -1))
   #RHOTHRES = -1.0E-6
   ϵ0 = 0.1E-12
   ϵ1 = 0.1E-16
   LIMIT = 50
   λlm = 0.0E+03
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
   H = zeros(Float64,length(perm),length(perm)) #Hessian
   #D = zero(H) #Elliptical Trust Region
   J = build_jcbn!(J,inds,vecs,params,perm)
   H, jtw = build_hess(jtw,J,W)
   uncs = paramunc!(uncs,H)
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
      βf,λlm,nomc,nrms,vals,nvecs,nparams = lbmq_turducken!(H,jtw,omc,λlm,nlist,inds,copy(params),perm,ofreqs)
      check = abs(nrms-rms)/rms
      if nrms < rms
         println(βf)
         rms = nrms
         omc .= nomc
         params .= nparams
#         println(params)
         #vecs .= nvecs
         J = build_jcbn!(J,inds,nvecs,params,perm)
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
         λlm /= 30.0
         if λlm ≤ 1.0E-24
            λlm = 0.0
         end
         stoit = 0
      else
         #params .= oparams
         λlm = max(2.0*λlm,1.0E-24)
      end
      #if rms < ϵ0
      #   converged = true
      #end
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
         break
      elseif (norm(βf))<ϵ1*(norm(params[perm])+ϵ1)
         slλ = (@sprintf("%0.4f", log10(λlm)))
         println("It would appear step size has converged. log₁₀(λ) = $slλ")
         break
      elseif counter ≥ LIMIT
         println("Alas, the iteration count has exceeded the limit")
         #println(omc)
         break
      else
         #write update to file
      end #check if
   end#while
   uncs = paramunc!(uncs,H)
   println(uncs)
   return params, vals
end
