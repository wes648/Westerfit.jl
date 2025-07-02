
function λgen(μ::Float64,er::Float64)::Float64
   ρ = 0.5
   er /= 1000.0 #convert err to GHz
   λ = ρ*μ*er #λF
   λ += (1.0 - ρ)*μ*er/(1+er) #λARC
   return λ
end
function lbmq_gain(β,λ::Float64,jtw,h,omc,nomc)::Float64
   out = 2β' * (λ*Diagonal(h)*β - jtw*omc)
   if out < 0
      printstyled("fucking gain function\n",color=:light_cyan)
   end
#   out = 2.0*(sum(abs2, omc .- nomc)) / out#abs(out)
   out = 2.0*(sum(abs2, omc) - sum(abs2, nomc)) / out#abs(out)
   return out
end
function lbmq_gain2(β,J,omc,nomc)::Float64
   pred = sum(abs2, J*β + omc)
   actu = sum(abs2,nomc)
   curr = sum(abs2,omc)
   #@show curr - pred > 0.0 
   if curr - pred < 0.0
      @warn "fucking gain function"
   end
#   @show curr - actu > 0.0 
   out = (curr - actu) / (curr - pred)
   return out
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
   Threads.@threads for i in 1:size(cfreqs,1)
      cfreqs[i] = vals[inds[i,3],inds[i,2]] - vals[inds[i,6],inds[i,5]]
   end
   #println(cfreqs)
   omc = ofreqs - cfreqs
   rms = BLAS.nrm2(omc)/√length(omc)
   #rms = norm(omc)/length(omc)
   return rms, omc, cfreqs
end
function rmscalc(cfreqs,ofreqs)
   omc = ofreqs - cfreqs
   rms = BLAS.nrm2(omc)/√length(omc)
   return rms, omc
end

function hstab(hssn,params)
   for i in 1:length(params)
      check = abs(hssn[i,i] / params[i])
      if check < 1e-12
         temp = abs(params[i])
         hssn[:,i] .*= temp
         hssn[i,:] .*= temp
      end
   end
end

function tsrapprox(j,β,stg)::Vector{Float64}
   out = zeros(size(j))
   out = j* β[collect(1:length(β))[(vcat(ones(18),stg) .> 0)]]
   return out#vec(sum(j,dims=2))
end

function paramunc(H,W,perm,omc)
   uncs = zeros(size(H,1))
   try
      uncs = diag(inv(Symmetric(H)))
   catch 
      uncs = diag(inv(H))
   end
   uncs .*= (omc' * W * omc)/(length(omc)-length(perm))
   return □rt.(uncs)
end
function correl(H)
   out = zeros(size(H))
   for i in 1:size(H,1), j in i:size(H,2)
      out[i,j] = H[i,j] / √(H[i,i]*H[j,j])
   end
   return Symmetric(out)
end
function covarr(corr,pσ)
   out = zeros(size(corr))
   for i in 1:size(corr,1), j in i:size(corr,2)
      out[i,j] = corr[i,j] * pσ[i] * pσ[j]
   end
   return Symmetric(out)
end
function covarr2(hess,omc)
   σ2 = 2*sum(abs2,omc)/(size(omc,1) - size(hess,1))^2
   return σ2 .* inv(hess)
end

function permdeterm(scls,stgs)
   out = collect(1:length(scls))[(scls .> 0) .* (vcat(ones(11),stgs) .≥ 0)]
end

function fincheck!(conv,endp,rms,βf,λlm,goal,check,ϵ0,ϵ1,ϵ2,counter,LIMIT,prms,grad)
   if (rms ≤ goal)#&&(counter > 1)
      println("A miracle has come to pass. The fit has converged")
      endp = "converge"
      conv = true
   elseif (check < ϵ0)
      println("The RMS has stopped decreasing. Hopefully it is low")
      endp = "RMS"
      conv = true
#   elseif (norm(βf))<ϵ1*(norm(prms)+ϵ1)
   elseif norm(βf ./ prms)<ϵ1
   #This stopping criteria needs to be scaled for the wildly varying parameter
   #magnitues. 3 dec 24
      slλ = (@sprintf("%0.4f", log10(λlm)))
      println("It would appear step size has converged. log₁₀(λ) = $slλ")
      @show norm(βf)
      endp = "step size"
      conv = true
   elseif (λlm > 1.0e+9)#&&(Δlm == 0.0)
      println("λlm exceeded threshold.")
      println("If you were using the turducken, try again without it")
      endp = "LMthresh"
      conv = true
   elseif norm(grad) < ϵ2
      println("Gradient is now quite small! This should be good")
      endp = "grad"
      conv = true
   elseif counter ≥ LIMIT
      println("Alas, the iteration count has exceeded the limit")
      endp = "iter"
      conv = true
   else
   end #check if
   return conv, endp
end
