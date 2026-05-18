"""
TODO:
 - first derivatives: 
 - Hessian JtWJ: 
 - Hessian Birss: 
 - Hessian Magnus: 


potential derive struct:
   Array{LowerTriangular{Float64,3},2}
   big array indecese are J & σ
   subarray are [nb, nk, i] were in nb is ⟨n'| , nk is |n⟩, and i is the ith operator
   the elements are ⟨n'| dℋ/di |n⟩
"""



#derivative can be calculated from :
# enact(O::Op,ψ::Psi,wvs::Eigs, 1.0, UR, true)

function sumder()
   ind = rpid+1
   if ind ≤ length(stg)+ HCLENGTH
      check = stg[ind-HCLENGTH]
      while check < zero(check)
         pm = prm[ind]
         out .+= enact(O[ind],ψ,wvs, prm[ind-HCLENGTH, UR, true)
         #if sum(ms .% 3)≈0 && j==1.5
         #   @show out
         #end
         ind += 1
         if ind-HCLENGTH ≤ length(stg)
            check = stg[ind-HCLENGTH]
         else
            check = 0
         end
      end
   end
   return out
end
function derivmat(wvs)
   if scl[rpid] < 0 #should this be ≤ 0 ???
   elseif rpid ≤ 4 #pure rot
   elseif 5 ≤ rpid ≤ 9 #spin-rot
   elseif 10 ≤ rpid ≤ 12 #qua
   elseif rpid==13 # F
   elseif rpid==16 # Vnf
   elseif rpid==14 # ρzF
   elseif rpid==15 # ρxF
   elseif rpid==17 # ηz
   elseif rpid==18 # ηx
   else #user def
      out = enact(O::Op,ψ::Psi,wvs::Eigs, 1.0, UR, true)
      out .= sumder()
   end
   return out
end

function anaderiv(...)
   mat = derivmat(...)
   out = sand(mat, wvs.rst.vecs[] )
   return diag(out)
end

function der2_magnus_elem(dx,dy,h,e,v)
   2*sand(dx*pinv(e*I(size(h,1)) - h)*dy, v)
end
function der2_birss_elem(dx,dy,wvs,i)
   e_i = wvs.rst.vals[lims, i] # <----------
   vec_i = wvs.rst.vecs[lims, i] # <----------
   inds = filter(x->!isequal(x,i), lims)
   out = 0.0
   for x ∈ inds
      vec_j = wvs.rst.vecs[lims, x] # <----------
      e_j = wvs.rst.vals[lims, j] # <----------
      out += (vec_i' * dx * vec_j * vec_j' * dy * vec_i) / (e_i - e_j)
   end
   return 2.0*out
end
