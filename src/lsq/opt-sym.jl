"""
This is a test document to better understand LBMQ optimization.
I am going to apply it to a symmetric top E(N,K) = (A-B)*K^2 + B*N*(N+1)
J(N,K) = [K^2; N*(N+1)-K^2]
Having two variables makes it plottable.
The true energy levels will be calculated 
  from E(N,K) = (A-B)*K^2 + B*N*(N+1) - DN*N^2*(N+1)^2 - DK*K^4 - DNK*N*(N+1)*K^2
This will provide error for additional testing
"""

using LinearAlgebra, Plots, Printf

ntop(n) = n*(n+1)

function Enk(p,q)
   out = zeros(size(q,1))
   @. out = p[2] * ntop(q[:,1])
   @. out += (p[1] - p[2])*q[:,2]^2
   return out
end
function Enk_cd(p,q)
   out = zeros(size(q,1))
   @. out = p[2] * ntop(q[:,1])
   @. out += (p[1] - p[2])*q[:,2]^2
   @. out += p[3] * ntop(q[:,1])^2
   @. out += p[4] * q[:,2]^4
   @. out += p[5] * ntop(q[:,1])*q[:,2]^2
   return out
end

function Jnk(p,q)
   out = zeros(size(q,1),2)
   @. out[:,1] = q[:,2]^2
   @. out[:,2] = ntop(q[:,1]) - q[:,2]^2
   return out
end
function build_hess(J,W)
   jtw = J' * W
   h = jtw * J
   return h, jtw
end

function qngen(n)
   out = zeros(Int,n+1,2)
   out[:,1] .= n
   out[:,2] = collect(0:n)
   return out
end
function qngen_tot(nmax)
   out = qngen(0)
   for n in 1:nmax
      out = vcat(out,qngen(n))
   end
   return out
end

function λgen(μ::Float64,er::Float64)::Float64
   ρ = 0.5
   er /= 1000.0 #convert err to GHz
   λ = ρ*μ*er #λF
   λ += (1.0 - ρ)*μ*er/(1+er) #λARC
   return λ
end
function lbmq_gain(β,λ::Float64,g,h,omc,nomc)::Float64
   out = transpose(β)*(λ*Diagonal(h)*β - g)
   if out < 0
      println("fucking gain function")
   end
   out = 2.0*(sum(omc .^2) - sum(nomc .^2)) / out#abs(out)
   return out
end
function fincheck!(converged,endpoint,rms,βf,λlm,goal,check,ϵ0,ϵ1,counter,LIMIT,prms)
   if (rms ≤ goal)#&&(counter > 1)
      println("A miracle has come to pass. The fit has converged")
      endpoint = "converge"
      converged = true
   elseif (check < ϵ0)
      println("The RMS has stopped decreasing. Hopefully it is low")
      endpoint = "RMS"
      converged = true
   elseif (norm(βf))<ϵ1*(norm(prms)+ϵ1)
      slλ = (@sprintf("%0.4f", log10(λlm)))
      println("It would appear step size has converged. log₁₀(λ) = $slλ")
      @show norm(βf)
      endpoint = "step size"
      converged = true
   elseif (λlm > 1.0e+9)#&&(Δlm == 0.0)
      println("λlm exceeded threshold.")
      println("If you were using the turducken, try again without it")
      endpoint = "LMthresh"
      converged = true
   elseif counter ≥ LIMIT
      println("Alas, the iteration count has exceeded the limit")
      endpoint = "iter"
      converged = true
   else
   end #check if
   return converged, endpoint
end

function rmscalc(test,real)
   omc = real - test
   rms = BLAS.nrm2(omc)/√(length(omc))
   return rms, omc
end

function lbmq_step(β,H,J,jtw,omc,λ)
   A = Hermitian(H + λ*Diagonal(H))
   #while isposdef(A)==false #this could be tidier
   #   λ = max(2.0*λ,1.0E-24)
   #   A = Hermitian(H + λ*Diagonal(H))
   #   printstyled("yike\n"; color=:green)
   #end
   @show A
   if isinf(λ)
      @warn "LB-MQ Matrix not pos-def!"
   end
   A = cholesky!(A)
   β .= ldiv!(β, A, jtw*omc)
   return β,λ
end

function lbmq_test()
   LIMIT = 10
   qns = qngen_tot(8)
   Etrue = Enk_cd([5.07677915; 0.84575606; 0.185784e-08; 0.649418e-07; 0.138846e-07], qns)
#   Etrue = Enk_cd([5.07677915; 0.84575606; 0.185784e-05; 0.649418e-04; 0.138846e-04], qns)
   #Etrue = Enk([5.07677915; 0.84575606], qns)

   params = [6.2366; 1.5143]
   #params = [5.02366; 0.8443]
   Etest = Enk(params,qns)
   uncs = fill(1.3e-7,length(Etest))

   rms, omc = rmscalc(Etest,Etrue)
   paramarray = zeros(2,LIMIT)
   paramarray[:,1] = params

   println("Initial RMS = $rms")
   goal = BLAS.nrm2(uncs)/√length(uncs)
   W = Diagonal(1.0 ./ uncs)
   ϵ0 = 0.0E-6 #RMS check
   ϵ1 = 0.1E-10 #step size check
   μlm = 0.0001
   λlm = λgen(μlm, rms) 
   println("Initial λ = $λlm")
   counter = 0

   nparams = copy(params)
   βf = zeros(length(params)) #step
   J = Jnk(params,qns)
   @show size(J)
   H, jtw = build_hess(J,W)
   @show H
   if true ∈ isnan.(H)
      println("FUCKING FUCKING FUCK. NaN in Hessian")
   end
   endpoint = "not yet"
   converged=false
while converged==false
   @show -2jtw*omc

   λlm = λgen(μlm, rms) 
   βf,λ = lbmq_step(βf,H,J,jtw,omc,λlm)
   nparams = params + βf

   Etest = Enk(nparams,qns)
   nrms, nomc = rmscalc(Etest,Etrue)
   ρlm = lbmq_gain(βf,λlm,-2jtw*omc,H,omc,nomc)
   @show βf
   @show ρlm
   @show nrms
   check = abs(nrms-rms)/rms

#   if (ρlm > 1.0e-7)#||(nrms<rms)#1.0e-7
   if nrms < rms
      rms = nrms
      params .= nparams
      if true#(mod(counter,2)==0)&&(counter > 0)
         #println("New vector time!")
         J = Jnk(params,qns)
         H, jtw = build_hess(J,W)
      else
         y = -2*jtw*(nomc - omc)
         H += (y*y')/(y' *βf) - (H*βf*βf' *H')/(βf' *H*βf)
#         cfreqs = tsrapprox(J,params)
      end
      omc = nomc
      counter += 1
      paramarray[:,counter+1] = params
      #@show norm(jtw*omc)
      μlm /= 30.0
   else
      printstyled("grumble grumble\n"; color=:red)
      μlm = max(4.0*μlm,1.0E-24)
   end
   converged,endpoint = fincheck!(converged,endpoint,rms,βf,λlm,goal,check,ϵ0,ϵ1,counter,LIMIT,params)
   if converged
      @show βf
      @show J
      @show -2jtw*omc
      @show eigvals(H)
   end
end#while
   #puncs = zeros(size(params))
   #puncs[perm] = paramunc(H,W,perm,omc)
   #covarmat = covarr(correl(H),puncs)
   @show paramarray
end
