
angle(a::Vector,b::Vector)::Float64 = dot(a,b) / (norm(a)*norm(b))

function λgen(μ::Float64,er::Float64)::Float64
   ρ = 0.5
   er /= 1000.0 #convert err to GHz
   λ = ρ*μ*er #λF
   λ += (1.0 - ρ)*μ*er/(1+er) #λARC
   return λ
end

function jlister(inds::Matrix{Int})::Matrix{Int}
   #finds all the unique J & σ pairs
   js = vcat(inds[:,1],inds[:,4])
   σs = vcat(inds[:,3] .+1,inds[:,6] .+1)
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

function χ2calc(wvs::Eigs,lins::Lines, lnjct)#::Float64
   cfreqs = zero(lins.frqs)
   Threads.@threads for i in 1:size(cfreqs,1)
      cfreqs[i] = wvs.rst.vals[lins.inds[i,2],lins.inds[i,3]+1] - 
                        wvs.rst.vals[lins.inds[i,5],lins.inds[i,6]+1]
   end
   omc = lins.frqs - cfreqs
   χ2 = sand(Diagonal(lins.wght[1:end .!= lnjct]).^2, omc[1:end .!= lnjct])
   return χ2, omc, cfreqs
end
function rmscalc(omc::Vector{Float64}, ℋ, prjct::Vector{Int}, lnjct::Vector{Int}, χ2::Float64)
   nparam = sum(x -> x.scl, ℋ) - length(prjct) + 1
   dof = length(omc) - length(lnjct) - nparam + 1
   rms = √(sum(abs2, omc)/dof)
   @show dof
   return rms, √(χ2/dof)
end

function correl(H,wrms)
   out = inv(H)
   uncs = □rt.(diag(out) .* wrms^2)
   for i in 1:size(H,1), j in i:size(H,2)
      out[i,j] = out[i,j] / (uncs[i]*uncs[j])
#      out[i,j] = out[i,j] / (√(out[i,i]*out[j,j]))
   end
   return Symmetric(out), uncs
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

function stepsizehcekcer(δ,ℋ,prjct)
   out = 0.0
   j = 1
   for i ∈ 1:length(ℋ)
      if !iszero(ℋ[i].scl) && j ∉ prjct
         out += (δ[j] / ℋ[i].val)^2
         j += 1
      end
   end
   return √out
end

function fincheck(ctrl,βf,λlm,check,counter,ℋ,prjct,grad,wrms)
   ϵ0 = 0.1E-8 #rms change threshold
   ϵ1 = 0.1E-6 #step size threshold
   ϵ2 = 0.1E-3 #gradient threshold
   if (wrms ≤ ctrl.goal)#&&(counter > 1)
      println("A miracle has come to pass. The fit has converged")
      endp = "converge"
      conv = true
   elseif (check) < ϵ0
      @show check
      println("The RMS has stopped decreasing. Hopefully it is low")
      endp = "RMS convergence"
      conv = true
#   elseif (norm(βf))<ϵ1*(norm(prms)+ϵ1)
   elseif stepsizehcekcer(βf,ℋ,prjct) <ϵ1
   #This stopping criteria needs to be scaled for the wildly varying parameter
   #magnitues. 3 dec 24
      slλ = (@sprintf("%0.4f", log10(λlm)))
      println("It would appear step size has converged. log₁₀(λ) = $slλ")
      @show norm(βf)
      endp = "step size convergence"
      conv = true
   elseif (λlm > 1.0e+9)#&&(Δlm == 0.0)
      println("λlm exceeded threshold.")
      println("If you were using the turducken, try again without it")
      endp = "λlm exceeded threshold"
      @show λlm
      conv = true
   elseif norm(grad) < ϵ2
      println("Gradient is now quite small! This should be good")
      endp = "gradient convergence"
      conv = true
   elseif counter ≥ ctrl.maxiter
      println("Alas, the iteration count has exceeded the limit")
      endp = "maximum number of iterations achieved"
      conv = true
   else
      conv = false
      endp = ""
   end #check if
   return conv, endp
end
