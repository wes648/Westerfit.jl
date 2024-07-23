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

function autsugdis()
end

function perm_asd(s,nfold)
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

function eder_asd(J)
   der = zeros(size(J))
   for i in 1:size(J,2)
      @views der[:,i] = log.(J[:,i]) .* J[:,i]
   end
   return der
end

function pred_asd(J,p)
   out = zeros(size(J,1))
   out .= prod(J .^ p', dims=2)
   return normalize!(out)
end