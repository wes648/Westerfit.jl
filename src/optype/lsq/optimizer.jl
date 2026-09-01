
function taylor_approx(δ::Vector{Float64},omc::Vector{Float64},
                       J::Matrix{Float64},H::Matrix{Float64},prjct,lnjct)::Float64
   Jt = J[1:length(omc) .!= lnjct, 1:size(H,1) .!= prjct]
   Ht = H[1:size(H,1) .!= prjct, 1:size(H,1) .!= prjct]
   omct = omc[1:length(omc) .!= lnjct]
   return (δt' * Jt' * omct) + 0.5*sand(Ht,δt)
end

function prm_rjct(H)::Vector{Int}
   # if hessian diagonal element is too small, temporarily remove from H & J
   perm = zeros(Int,1)
   dH = diag(H)
   @inbounds for i ∈ eachindex(dH)
      if abs(dH[i]) < TOL
         perm = vcat(perm,i)
      end # if
   end # for
   #if isempty(perm); perm = [0]; end
   return perm
end

function lin_rjct(J,W,omc)::Vector{Int}
   # if omc/W > THRESH, temporarily remove from J & W & omc
   perm = zeros(Int,1)
   check = (omc ./ W).^2
   @inbounds for i ∈ eachindex(check)
      if check[i] < TOL
         perm = vcat(perm,i)
      end
   end
   #if isempty(perm); perm = [0]; end
   return perm
end

function hess_prep(H,λ,prjct)
   A = H[1:end .!=prjct, 1:end .!=prjct]
   A = Hermitian(A + λ*I(size(A,1)))
   while !isposdef(A) #this could be tidier
      λ = max(2.0*λ,1.0E-24)
      A += (λ - 0.5*λ)*I(size(A,1))
      if isinf(λ)
      @warn "Hessian Matrix not pos-def!
      Make sure you aren't trying to optimize a parameter with value of 0.0.
      This code is about to crash"
         break
      end
   end
   return cholesky!(A), λ
end
function grad_prep(J,W,omc)
   perm = lin_rjct(J,W,omc)
   rjct = 1:length(omc) .!= perm
   grad = -J[rjct,:]' * W[rjct] * omc[rjct]
   return grad, perm
end

function newt_step(β,λ,Δ,H,J,w,omc)
   prjt = prm_rjct(H)
   grad, lrjt = grad_prep(J[:, 1:end .!=prjt], w, omc)
   A, λ = hess_prep(H,λ, prjct)
   β = ldiv!(A, grad)
   β = trust_region(β,A,Δ)
   return β,λ, prjct, lrjt
end

function update_prm(ℋ::Vector{Op},δ::Vector{Float64},prjct::Vector{Int})::Vector{Op}
   j = 1
   @inbounds for i ∈ eachindex(ℋ)
      if !iszero(ℋ[i].scl) && j ∉ prjct
         ℋ[i].val += δ[j]
      end
      j += 1
   end
   return ℋ
end
function trust_region(δ,H,Δ)
   if norm(B*δ) > Δ
      δ *= 0.5*Δ/norm(B*δ)
   end
   return δ
end


function opt_calc(molnam::String,ctrl::Controls,ℋ::Vector{Op},lins::Lines)
   output_init(molnam::String,ctrl::Controls,ℋ::Vector{Term})
   converged = false
   nwvs = Eigs(ctrl)
   wvs = Eigs(ctrl)
   jsσs = jlister(lins.inds)
   H_calc(ctrl,wvs,ℋ,jsσs)
   χ2, omc, cfq = χ2calc(wvs,lins)
   lχ2 = copy(χ2)
   J,H = deriv_calc(ctrl,ℋ,wvs,lins, γ)
   δ = zeros(nfit_count(ℋ)); oδ = ones(nfit_count(ℋ))
   println("Intitial wrms = ", √(χ2/length(omc)))
   counter = 0

   Δ = ctrl.Δlm0
   μlm = ctrl.λlm0

   while !converged
      λ = λgen(μlm, χ2)
      δ,λ,prjct,lrjt = newt_step(δ,λ,H,J,lins.W,omc)
      θ = round(angle(δ,oδ),3)
      approx = χ2 + taylor_approx(δ,omc,J,H, prjct,lrjt)

      ℋ = update_prm(ℋ,δ)
      H_calc(ctrl,nwvs,ℋ,jsσs)
      nχ2, nomc, ncfq = χ2calc(nwvs,lins)
      ρdn = (χ2 - nχ2) / (χ2 - approx)
      if (BOLD == 0)|| χ2 > 1e5
         stepcheck = (( nχ2 < χ2 )&&( ρdn >1e-6))
      else
         stepcheck = (( nχ2 < χ2 )&&( ρdn >1e-6)) || ((nχ2*(1-θ)^BOLD)<0.2*lχ2)
      end
     
      if stepcheck
         lχ2 = min(nχ2, lχ2)
         wvs = copy(nwvs)
         oδ = copy(δ)
      @time deriv_calc!(J,H, ctrl,ℋ,wvs,lins, γ)
         χ2, omc, cfq = χ2calc(wvs,lins)
         printstyled(" accepted step: ", color=:green)
         println("χ2 = ",lpad(χ2,5), ", ||δ|| = ", lpad(norm(δ),3), ", θ = " lpad(θ,3))
         sum_print(io, ℋ, prjct, lnjct, rms, wrms, χ2)
         μlm /= 10.0
         Δlm *= 2.0
      else
         # model failed to decrease
         printstyled(" rejected step: ", color=:red)
         println("χ2 = ",lpad(χ2,5), ", ||δ|| = ", lpad(norm(δ),3), ", θ = " lpad(θ,3))
         μlm *= 4.0
         Δlm *= Δlm
         ℋ = update_prm(ℋ,-δ)

      end # stepcheck if
      converged, enp = fincheck(ctrl,δ,λlm,nχ2-χ2,counter,ℋ[:].val, J'*W*omc)
   end #converged loop
   opt_finalize(molnam)
   return ℋ, wvs
end #function

function opt_finalize(molnam)
   corr = h2conv(H)
   output_final()
end




