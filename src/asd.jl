"""
This is the Automated Suggestor of Distortions
   It works by calculating the expectation values of the allowed 2nd order 
   operators & optimizing their powers. The energy contribution is approxmiated
   as ⟨A⟩ᵃ⟨B⟩ᵇ⟨C⟩ᶜ where A, B, and C are operators. Then the differences are
   calculated for each energy level and the normalized contributions are 
   compared to the normalized omcs. The normamlization is done to not have to
   also optimize a parameter. The %1.2f operators and the norm(omc)/norm(cont)
   will be returned in the format of a westerfit input line for the user to 
   decide how to implement.
"""

"""
perm_asd(s::Float64,nfold::Int64)
operator permuation generator for ASD
This determines which second order operators are used in the energy approximator
"""
function operm_asd(s,nfold)
   if s==zero(s) && nfold == zero(nfold)
      perm = collect(1:3)
   elseif s==0.5 && nfold == zero(nfold)
      perm = collect(1:9)
   elseif s≥1.0 && nfold == zero(nfold)
      perm = collect(1:12)
   elseif s==zero(s) && nfold ≠ zero(nfold)
      perm = [collect(1:4); collect(13:16)]
   elseif s==0.5 && nfold ≠ zero(nfold)
      perm = [collect(1:9); collect(13:18)]
   elseif s≥1.0 && nfold ≠ zero(nfold)
      perm = collect(1:18)
   else
      println("This case shouldn't occurr but I'll just turn on everything")
      perm = collect(1:18)
   end
   return perm
end

function jcb_asd(jlist,inds,ctrl,vecs,perm,scals,stg,tvecs)
   lst = vcat(inds[:,1:3],inds[:,4:6])
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   jcbn = zeros(Float64,size(inds,1),length(perm))
   deriv = derivcalc_2stg(jlist,ctrl,perm,vecs,tvecs)
   @threads for p in 1:length(perm)
   @simd for a in 1:size(inds,1)
      jcbn[a,p] = deriv[lst[a,3], lst[a,2]+1,p]
   end
   end
   return jcbn
end

"""
d/dx a^x = log(a)*a^x
energy derivative for ASD
"""
function eder_asd(J,p)
   #der = zeros(size(J))
   der = log.(J) .* (abs.(J).^p') .* sign.(J).^(isodd.(round.(p')))
   return der
end

"""
Prediction calculator for ASD
We are using a^x ≈ 0.5 |a|^x * cos(π x) (1 - sgn(a)) 
"""
function pred_asd(J,p)
   l = floor(Int,0.5*size(J,1))
   out = abs.(J) .^ p' .* cosp.(p) .* (0.5 .- 0.5.*sign.(J))
   out = prod(out, dims=2)
   return normalize!(out[1:l] .- out[l+1:end])
end
function func_asd(J,p,omc)
   pred = pred_asd(J,p)
   a = rmscalc( pred,omc)^2
   b = rmscalc(-pred,omc)^2
   if a < b
      return pred, 1
   else# b < a
      return -pred, -1
   end
end

function grad_asd(p,J,pred,omc,k)
   out = -2*k*(omc - pred)*eder_asd(J,p)
   return sum(out,dims=2)
end

normto(v::Array{Float64,1},x::Number) = x .* normalize!(v)

function autsugdis()
   ps = zeros(length(perm),maxiter+1

end
