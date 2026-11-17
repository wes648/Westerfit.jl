
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
   nparam = sum(x -> !iszero(x.scl), ℋ) - length(prjct) + 1
   dof = length(omc) - length(lnjct) - nparam + 1
   if 0 > dof
      @warn "Negative degrees of freedom. time to die"
   end
   rms = √(sum(abs2, omc)/dof)
   return rms, √(χ2/dof)
end

function correl(H,wrms)
   out = inv(H)
   uncs = □rt.(diag(out) .* wrms^2)
   for i in 1:size(H,1), j in i:size(H,2)
#      out[i,j] = out[i,j] / (uncs[i]*uncs[j])
      out[i,j] = H[i,j] / (√(H[i,i]*H[j,j]))
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

function stepsizechecker(δ,ℋ,prjct)
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

function accept_print(counter,wrms,rms,δ,ℋ,prjct,θ,ρ)
   printstyled(" accepted step: ", lpad(counter,4), color=:green)
   println(" wrms = ", lpad(round(wrms, digits=5),5), ", rms = ", lpad(round(rms,digits=5),3))         
   println("ρ = ", lpad(round(ρ,digits=5),3),
      ", ‖δ‖ = ", lpad(round(stepsizechecker(δ,ℋ,prjct), digits=5),3),
      ", θ = ", lpad(round(θ,digits=5),3))
end
function reject_print(counter,wrms,rms,δ,ℋ,prjct,θ,ρ)
   printstyled(" rejected step: ", lpad("",4), color=:red)
   println("wrms = ", lpad(round(wrms, digits=5),5), ", rms = ", lpad(round(rms,digits=5),3))
   println("ρ = ", lpad(round(ρ,digits=5),3),#"χ2 = ",lpad(round(χ2, digits=5),5),", "
      ", ‖δ‖ = ", lpad(round(stepsizechecker(δ,ℋ,prjct), digits=5),3), 
      ", θ = ", lpad(round(θ,digits=5),3))
end


function fincheck(ctrl,βf,λlm,check,counter,ℋ,prjct,grad,wrms)
   ϵ0 = 0.1E-8 #rms change threshold
   ϵ1 = 0.1E-5 #step size threshold
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
   elseif stepsizechecker(βf,ℋ,prjct) <ϵ1
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

function conv_tr(ϵs::Vector{Float64},iters::Int,wrms::Float64,check::Float64, βf,ℋ,prjct,
      grad,Δ::Float64)#::Bool
   count = zeros(Int,length(ϵs)-1)
   if iters ≥ ϵs[2]
      println("Alas, the iteration count has exceeded the limit")
      count[1] = 1 + length(ϵ)
   end
   if wrms ≤ ϵs[3]
      println("A miracle has come to pass. The fit has converged")
      count[2] = 1
   end
   if check ≤ ϵs[4]
      @show check
      println("The RMS has stopped decreasing. Hopefully it is low")
      count[3] = 1
   end
   if stepsizechecker(βf,ℋ,prjct) ≤ ϵs[5] # <----- not exactly the write check function
      @show norm(βf)
      println("It would appear step size has converged.")#" log₁₀(λ) = $slλ")
      count[4] = 1
   end
   if norm(grad, Inf) ≤ ϵs[6]
      @show norm(grad, Inf)
      println("Gradient is now quite small! This should be good")
      count[5] = 1
   end
   if Δ ≤ ϵs[7]
      @show Δ
      println("Alas, the trust region has become too restrictive")
      count[6] = 1 + length(ϵ)
   end
   return (sum(count) ≥ ϵs[1] ? true : false), count
end

function step_accepter(ctrl, ρ, wrms, χ2, nχ2, lχ2)
   if (ctrl.BOLD == 0) || wrms > 5e2
      stepcheck = (( nχ2 < χ2 )&&( ρ > 1e-4))
   else
      stepcheck = (( nχ2 < χ2 )&&( ρ > ctrl.ρfit)) || nχ2*(1-θ)^ctrl.BOLD < 0.2*lχ2
   end
   return stepcheck && !iszero(ρ)
end
