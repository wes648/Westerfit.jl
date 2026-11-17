
# Dogleg implementation based on dogleg.rs

normsq(x)::Float64 = sum(abs2, x)

function τdl_calc(Δ, psd, pgn)::Float64
   ξ = δsd' * (δgn - δsd) / normsq(δgn - δsd)
   τ = 1 - ξ + √((Δ^2 - normsq(δsd))/normsq(δgn - δsd) + ξ^2)
   return τ
end

function damp_init(J)::Diagonal{Float64,Vector{Float64}}
   out = ones(size(J,2))
   for i ∈ eachindex(out)
      out[i] += norm(J[:,i])
   end
   return Diagonal(out)
end
function damp_update(D,J)::Diagonal{Float64,Vector{Float64}}
   out = zero(D)
   for i ∈ 1:size(D,1)
      out[i,i] = max(min(1e16, norm(J[:,i]) ), 1e-3)
   end # for 
   return D*out
end

function grad_step(J,g)
   #sign here might not be correct shoudld chekc
   return normsq(g)/normsq(J*g) .* g
end
function hess_prep(H,D,λ)
   #μ = 1e-8
   A = Hermitian( H + λ*D^2 )
   while !isposdef(A)
      μ *= 10.0
      A = Hermitian( H + λ*D^2 )
      if μ > 100.0
         @warn "Regularization Parameter exceeded threshold. Hessian is not posdef!
         Program will now crash"
            break
      end
   end
   return cholesky!(A)
end
function hess_step(H,D,grad,λ)
   A = hess_prep(H,D,λ)
   return ldiv!(A, grad)
end

function dogleg_scaled(Ds, J, grad, Hess, Δ,λ)::Vector{Float64}
   D = damp_update(Ds,J)
   δsd = grad_step(J,grad)
   if norm(D*δsd) ≥ Δ
      return (Δ/norm(δsd)) .* δsd
   else
      δgn = hess_step(Hess, D, grad,λ)
      if norm(D*δgn) ≤ Δ
         return δgn
      else
         τ = τdl_calc(Δ, D*δsd, D*δgn)
         if τ ≤ 1.0
            δdl = τ*δsd
         else # τ ∈ (1,2]
            δdl = δsd + (τ-1)*(δgn - δsd)
         end # τ if
         return δdl
      end # gn if
   end #sd if
end

function Δ_update(Δ,p,ρ)
   if ρ > 0.75
      Δ = max(Δ, 3*norm(p) )
   elseif ρ < 0.25
      Δ *= 0.5
   else
      #nothing
   end
   return Δ
end

function gain_dl(χ2,nχ2, γ, grad, p, H)
   denom = (p'*grad + 0.5*sand(H,p))
   if denom < 0.0
      @warn "FUCK GAIN DENOM NOT POSITIVE"
   end
   numer = (χ2 - nχ2)
#   @show denom
#   @show numer
   return numer / denom
end


function μ2λlm(μ,wrms)
   if iszero(μ)
      return 1e-8
   else
      return wrms > 1e3 ? μ : λgen(μ, wrms)
   end
end

function opt_calc_dl(molnam::String,ctrl::Controls,ℋ::Vector{Term},lins::Lines)
   output_init(molnam,ctrl,ℋ,lins)
   nwvs = Eigs(ctrl);   wvs = Eigs(ctrl)
   jsσs = jlister(lins.inds)
   @time "energy time" H_calc(ctrl,wvs,ℋ,jsσs)
   χ2, γ, cfq = χ2calc(wvs,lins, [0]); lχ2 = 2*χ2
   @time "deriv time" J,H = deriv_calc(ctrl,ℋ,wvs,lins, γ)
   Ds = damp_init(J)

   δ = zeros(nfit_count(ℋ)); oδ = ones(nfit_count(ℋ))
   sum_init_print(molnam, sum(map(x->x.scl, ℋ)), lins, γ)
   prjct = [0]; lrjt = [0]
   rms,wrms = rmscalc(γ, ℋ, prjct, lrjt, χ2)
   println("Intitial wrms = ", wrms)

   counter = 0
   Δ = ctrl.Δlm0;   μlm = ctrl.λlm0
   ϵs = [1; ctrl.maxiter; ctrl.goal; 1e-9; 1e-6; 1e-4; 1e-12] #convergence criteria
   # min count, iter, goal, rms change, stepsize, grad, Δ
   converged = false
   how = zeros(Int,length(ϵs)-1)

   while !converged
      @show Δ
      λ = μ2λlm(μlm,wrms)
      δ = dogleg_scaled(Ds, J, J'*γ, H, Δ,λ)
      θ = !iszero(counter) ? min(round(angle(δ,oδ),digits=3), 0.9999) : 0.0

      ℋ = update_prm(ℋ,δ,prjct)
      @time "energy time" nwvs = H_calc(ctrl,nwvs,ℋ,jsσs)
      nχ2, nγ, ncfq = χ2calc(nwvs,lins,lrjt)
      nrms,nwrms = rmscalc(γ, ℋ, prjct, lrjt, nχ2)
      check = □rt(abs(nχ2 - χ2)/χ2)

      ρ = gain_dl(χ2, nχ2, γ, J'*γ, δ, H)
#      @show ρ
#      @show nχ2
#      @show χ2
      Δ = Δ_update(Δ,δ,ρ)

      if step_accepter(ctrl, ρ, wrms, χ2, nχ2, lχ2)
         # good step!
         counter += 1
         lχ2 = min(nχ2, lχ2)
         wvs = deepcopy(nwvs)
         oδ .= δ
      @time "deriv time" deriv_calc!(J,H, ctrl,ℋ,wvs,lins, γ)
         χ2, γ, cfq = χ2calc(wvs,lins,lrjt)
         rms,wrms = rmscalc(γ, ℋ, prjct, lrjt, χ2)
         accept_print(counter,wrms,rms,δ,ℋ,prjct,θ,ρ)
         sum_print(molnam, counter, ℋ,δ, prjct, lrjt, rms, wrms, χ2)
         μlm /= 20.0
      else
         # bad step!
         reject_print(counter,nwrms,nrms,δ,ℋ,prjct,θ,ρ)
         μlm *= 10.0
         ℋ = update_prm(ℋ,-δ,prjct)
      end # stepcheck if
      converged, how = conv_tr(ϵs, counter,wrms,check,δ,ℋ,prjct, J'*γ,Δ)
   end # converged loop
   opt_finalize2(molnam,ℋ,H,prjct,lrjt, rms, wrms, χ2, how)
   return ℋ, γ, χ2, cfq
end

function opt_finalize2(molnam,ℋ,H,prjct, lnjct, rms, wrms, χ2, how::Vector{Int})
   corr, unc = correl(H,wrms)
   output_final2(molnam,ℋ,unc,corr, prjct, lnjct, rms, wrms, χ2, how)
end