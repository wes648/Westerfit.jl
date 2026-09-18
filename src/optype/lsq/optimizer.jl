
function taylor_approx(δ::Vector{Float64},omc::Vector{Float64},
                       J::Matrix{Float64},H::Matrix{Float64},prjct,lnjct)::Float64
   Jt = @view J[1:length(omc) .!= lnjct, 1:size(H,1) .!= prjct]
   Ht = @view H[1:size(H,1) .!= prjct, 1:size(H,1) .!= prjct]
   omct = @view omc[1:length(omc) .!= lnjct]
   return (δ' * Jt' * omct) + 0.5*sand(Ht,δ)
end

function prm_rjct(H)::Vector{Int}
   # if hessian diagonal element is too small, temporarily remove from H & J
#   perm = zeros(Int,1)
#   dH = diag(H)
#   @inbounds for i ∈ eachindex(dH)
#      if abs(dH[i]) < 1e-3 # TOL
#         @warn "hey julia this Hessian tol should be controllable value"
#         perm = vcat(perm,i)
#      end # if
#   end # for
#   return perm
   return [0]
end

function lin_rjct(J,W,omc)::Vector{Int}
   # if omc/W > THRESH, temporarily remove from J & W & omc
#   perm = zeros(Int,1)
#   check = (omc ./ W).^2
#   @inbounds for i ∈ eachindex(check)
#      if check[i] > 1e5
#         @warn "hey julia this line tol should be controllable value"
#         perm = vcat(perm,i)
#      end
#   end
#   return perm
   return [0]
end

function hess_prep(H,λ,prjct)
   A = @view H[1:end .!=prjct, 1:end .!=prjct]
   A = Hermitian(A + λ*Diagonal(A) )
   while !isposdef(A) #this could be tidier
      λ = max(2.0*λ,1.0E-24)
      A += (λ - 0.5*λ)*Diagonal(A)
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
   grad = J[rjct,:]' * Diagonal(W[rjct]) * omc[rjct]
   return grad, perm
end

function newt_step(β,λ,Δ,H,J,w,omc)
   prjt = prm_rjct(H)
   grad, lrjt = grad_prep(J[:, 1:end .!=prjt], w, omc)
   A, λ = hess_prep(H,λ, prjt)
   β = ldiv!(A, grad)
   return β,λ, prjt, lrjt
end

function update_prm(ℋ::Vector{Term},δ::Vector{Float64},prjct::Vector{Int})::Vector{Term}
   j = 1
   @inbounds for i ∈ eachindex(ℋ)
      if !iszero(ℋ[i].scl) && j ∉ prjct
         ℋ[i].val += δ[j]
         j += 1
      end
   end
   return ℋ
end
function trust_region(δ,ℋ,prjct,Δ)
   #x = stepsizehcekcer(δ,ℋ,prjct)
   #@show x
   if norm(δ) > Δ #if norm(H*δ) > Δ
      δ *= Δ/norm(δ)
   end
   return δ
end


function opt_calc(molnam::String,ctrl::Controls,ℋ::Vector{Term},lins::Lines)
   output_init(molnam,ctrl,ℋ,lins)
   converged = false
   nwvs = Eigs(ctrl);   wvs = Eigs(ctrl)
   jsσs = jlister(lins.inds)
   @time "energy time" H_calc(ctrl,wvs,ℋ,jsσs)
   χ2, omc, cfq = χ2calc(wvs,lins, [0])
   endp = "not yet"
   lχ2 = copy(χ2)
   @time "deriv time" J,H = deriv_calc(ctrl,ℋ,wvs,lins, omc)
   δ = zeros(nfit_count(ℋ)); oδ = ones(nfit_count(ℋ))
   io = open(molnam*".out","a")
   sum_init_print(io, sum(map(x->x.scl, ℋ)), lins, omc)
   close(io)
   prjct = [0]; lrjt = [0]
   rms,wrms = rmscalc(omc, ℋ, prjct, lrjt, χ2)
   println("Intitial wrms = ", wrms)
   counter = 0
   Δ = ctrl.Δlm0 # 3e3 #
   μlm = ctrl.λlm0 #1e-6 #
   #λ = ctrl.λlm0
   while !converged

#      λ = wrms > 1e3 ? μlm : λgen(μlm, wrms)
      λ = μlm
      δ,λ,prjct,lrjt = newt_step(δ,λ,Δ,H,J, lins.wght, omc)
      δ = trust_region(δ,ℋ,prjct,Δ)
      θ = round(angle(δ,oδ),digits=3)
      approx = taylor_approx(δ,omc,J,H, prjct,lrjt)

      ℋ = update_prm(ℋ,δ,prjct)
      # @show δ
      @time "energy time" nwvs = H_calc(ctrl,nwvs,ℋ,jsσs)
      nχ2, nomc, ncfq = χ2calc(nwvs,lins,lrjt)
      ρdn = (χ2 - nχ2) / approx #(χ2 - approx)
      check = □rt(abs(nχ2 - χ2)/χ2)

      if (ctrl.BOLD == 0)|| √(χ2) > 1e5
         stepcheck = (( nχ2 < χ2 )&&( ρdn >1e-6))
      else
         stepcheck = (( nχ2 < χ2 )&&( ρdn >1e-6)) || ((nχ2*(1-θ)^ctrl.BOLD)<0.2*lχ2)
      end
     
      if stepcheck && !iszero(ρdn)
         counter += 1
         lχ2 = min(nχ2, lχ2)
         wvs = deepcopy(nwvs)
         oδ .= δ
      @time "deriv time" deriv_calc!(J,H, ctrl,ℋ,wvs,lins, omc)
         χ2, omc, cfq = χ2calc(wvs,lins,lrjt)
         rms,wrms = rmscalc(omc, ℋ, prjct, lrjt, χ2)
         printstyled(" accepted step: ", counter, color=:green)
         println(" wrms = ", lpad(round(wrms, digits=5),5), ", rms = ", lpad(round(rms,digits=5),3))         
         println(#"χ2 = ",lpad(round(χ2, digits=5),5),", "
            "||δ|| = ", 
            lpad(round(stepsizehcekcer(δ,ℋ,prjct), digits=5),3), ", θ = ", lpad(round(θ,digits=5),3))
         io = open(molnam*".out","a")
#         rms,wrms = rmscalc(omc, ℋ, prjct, lrjt)
         sum_print(io, counter, ℋ,δ, prjct, lrjt, rms, wrms, χ2)
         close(io)
         μlm /= 20.0
#         λ /= 20.0
         Δ *= 1.5
      else
         # model failed to decrease
         rms,wrms = rmscalc(omc, ℋ, prjct, lrjt, χ2)
         printstyled(" rejected step: ", color=:red)
         println("wrms = ", lpad(round(wrms, digits=5),5), ", rms = ", lpad(round(rms,digits=5),3))
         println(#"χ2 = ",lpad(round(χ2, digits=5),5),", "
            "||δ|| = ", 
            lpad(round(stepsizehcekcer(δ,ℋ,prjct), digits=5),3), ", θ = ", lpad(round(θ,digits=5),3))
         μlm *= 10.0
#         λ *= 10.0
         Δ *= 0.9
         ℋ = update_prm(ℋ,-δ,prjct)

      end # stepcheck if
converged, endp = fincheck(ctrl,δ,λ,check,counter, ℋ,prjct, J'*Diagonal(lins.wght)*omc,wrms)
   end #converged loop
   opt_finalize(molnam,ℋ,H,prjct,lrjt, rms, wrms, χ2, endp)
   close(io)
   return ℋ, omc, cfq
end #function

function opt_finalize(molnam,ℋ,H,prjct, lnjct, rms, wrms, χ2, endp)
   corr, unc = correl(H,wrms)
   output_final(molnam,ℋ,unc,corr, prjct, lnjct, rms, wrms, χ2, endp)
end




