
function build_jcbn!(jcbn,grad,inds,vecs,params,perm,omc)
"""
This builds the Jacobian based on the Hellmann–Feynman theorem.
"""
   Threads.@threads for a in 1:size(inds)[1]
      ju = 0.5*inds[a,1]
      jl = 0.5*inds[a,4]
      σu = inds[a,2]
      σl = inds[a,5]
      vecu = vecs[1:Int((2*S+1)*(2*ju+1)*(2*mcalc+1)),inds[a,3],σu+1]
      vecl = vecs[1:Int((2*S+1)*(2*jl+1)*(2*mcalc+1)),inds[a,6],σl+1]
      for i in 1:length(perm)
         b = perm[i]
         jcbn[a,i] = anaderiv(jl,S,σl,vecl,params,b) - anaderiv(ju,S,σu,vecu,params,b)
      end
   end
   grad = transpose(jcbn)*omc
   return jcbn, grad
end

function build_hess!(hssn,dk,jcbn,weights)
   hssn = transpose(jcbn)*weights*jcbn
   Threads.@threads for i in 1:size(hssn)[1]
      dk[i,i] = norm(hssn[:,i])
   end
   return hssn, dk
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
   β2 = zeros(Float64,size(β))
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
   β = ldiv!(β, A, -grad)
   return β,λ
end


function paramunc!(uncs,H)
   uncs = diag(inv(H))
   return uncs
end
function ellipictr()
end

function lbmq_opttr(nlist,ofreqs,uncs,inds,params,scales,λ)
   vals,vecs = limeigcalc(nlist, inds, params)
   rms, omc = rmscalc(vals, inds, ofreqs)
   perm,n = findnz(sparse(scales))
   println("Initial RMS = $rms")
   goal = sum(uncs)/length(uncs)*0.001
   newparams = copy(params)
   W = diagm(0=>(uncs .^ -1))
   converged = false
   THRESHOLD = 1.0E-8
   RHOTHRES = -1.0E-6
   ϵ0 = 0.1E-8
   ϵ1 = 0.1E-24
   LIMIT = 500
   λlm = 1.0E+03
   Δlm = 1.0E+1
   Δlm *= length(perm)
   counter = 0
   stoit = 5
   rms, omc = rmscalc(vals,inds,ofreqs)
   nparams = copy(params)
   uncs = zeros(Float64,size(perm))
   β = zeros(Float64,size(perm)) #step
   J = zeros(Float64,size(inds)[1],length(perm)) #Jacobian
   g = zeros(Float64,size(omc)) #gradient
   K = zeros(Float64,size(omc)) #acceleration correction
   #A = zeros(Float64,size(J)) #
   H = zeros(Float64,length(perm),length(perm)) #Hessian
   D = zeros(Float64,size(H)) #Elliptical Trust Region
   J, g = build_jcbn!(J,g,inds,vecs,params,perm,omc)
   H, D = build_hess!(H,D,J,W)
   uncs = paramunc!(uncs,H)
   converged=false
   while converged==false
      β,λlm = lbmq_step!(β,H,g,λlm)
      if false
         β = lbmq_acc!(K,β,J,W,nlist,inds,params,perm,omc,ofreqs,λlm)
      end
      β0 = β
      if (norm((H^(-1/2))*β)>Δlm)&&(λ!=0.0)
#      if (norm(β ./ params[perm])>Δlm)&&(λ!=0.0)
         β *= Δlm/norm(D*β)
      end
      nparams[perm] = params[perm] + β #+ (inv(H)^2)*β
      vals, nvecs = limeigcalc(nlist, inds, nparams)
      nrms, nomc = rmscalc(vals,inds,ofreqs)
      check = abs(nrms-rms)/rms
      if nrms < rms
         rms = nrms
         omc = nomc
         params = nparams
         vecs = nvecs
         J, g = build_jcbn!(J,g,inds,vecs,params,perm,omc)
         H, D = build_hess!(H,D,J,W)
         counter += 1
         srms = (@sprintf("%0.4f", rms))
         slλ = (@sprintf("%0.4f", log10(λlm)))
         sΔ = (@sprintf("%0.6f", Δlm))
         scounter = lpad(counter,3)
         println("After $scounter interations, RMS = $srms, log₁₀(λ) = $slλ, Δₖ = $sΔ")
         #println(β)
         #println(H^(-1/2))
         #println(β)
         #println(params[perm])
         λlm /= 3.0
         if λlm ≤ 1.0E-24
            λlm = 0.0
         end
         stoit = 0
      else
         λlm *= 2.0
         λlm = max(λlm,1.0E-24)
      end
      if rms < ϵ0
         converged = true
      end
      ρlm = lbmq_gain(β,λlm,g,rms,nrms)
      if ρlm ≥ 0.75
         Δlm *= 2.0
      else
         Δlm = max(0.90*Δlm,0.0001)
      end
      if (rms ≤ goal)#&&(counter > 1)
         println("A miracle has come to pass. The fit has converged")
         break
      elseif (check < THRESHOLD)
         println("The RMS has stopped decreasing. Hopefully it is low")
         break
      elseif (norm(β))<ϵ1*(norm(params[perm])+ϵ1)
         #if stoit ≤ 5
         #   println("Stocast!")
         #   λlm = 0.0
         #   params[perm] += (0.5 .- rand(length(perm))) .* uncs .* 4.0E+2
         #   stoit += 1
         #else
         println("It would appear step size has converged")
         break
         #end
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
