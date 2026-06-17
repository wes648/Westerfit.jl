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

Birss second derivatives:
S_xyn = d²E/dxdy @ state n = 2 ∑_n≠n' ⟨n| dℋ/dx |n'⟩⟨n'| dℋ/dy |n⟩ / E_n - E_n'
Taylor series: F(x) ≈ f(a) + df(a)/dx (x-a) + d²f(a)/dx² 1/2 (x-a)²
               F = f + (x-a)ᵀDf(a) + 1/2 (x-a)ᵀD²f(a)(x-a)


r_n = W/2(k - ν_ij)² =  W/2(k - E_i + E_j)²
χ² = γᵀ W γ

d r_n / dx = W(-dE_i/dx + dE_j/dx)*(k - E_i + E_j)
J_n = dE_i/dx - dE_j/dx
∇χ² = -Jᵀ W γ

d² r_n / dx² = W(-d²E_i/dx² + d²E_j/dx² )(k - E_i + E_j) + W(-dE_i/dx + dE_j/dx)²
d² r_n / dxdy = (-d²E_i/dxdy + d²E_j/dxdy )W(k - E_i + E_j) + W(-dE_i/dx + dE_j/dx)(-dE_i/dy + dE_j/dy)
S_xyn = d²E_i/dxdy - d²E_j/dxdy 
   = (2 ∑_k≠i ⟨i| dℋ/dx |k⟩⟨k| dℋ/dy |i⟩ / E_i - E_k') - (2 ∑_k≠j ⟨j| dℋ/dx |k⟩⟨k| dℋ/dy |j⟩ / E_j - E_k'))
H_xy = ∑_n d² r_n / dxdy 
H(χ²) = -∑_n S_xyn W_n γ_n + J'WJ

geodesic accel:
solve J α = -fvv
fvv is a vector of v' (S[:,:,n]) v
S[:,:,n] is 2nd deriv mat for transition n
so v' S[:,:,n] v is a number and [v' S[:,:,n] v] is a vector
δ = v + 0.5 α

init δ = zeros(length(perm))
init γ = zeros(length(ofreqs))
init W = Diagonal( unc .^-2 )
init J = zeros(length(ofreqs), length(perm))
init S = zeros(length(perm), length(perm), length(ofreqs))
init H = zeros(length(perm), length(perm))
update γ = ofreqs - cfreqs
update J = dE_i/dx - dE_j/dx
update S = d²E_i/dxdy - d²E_j/dxdy
update H = -sum(x->S[:,:,x]*W[x,x]*γ[x],eachindex(γ)) + J'WJ
update δ = (H + λ*Diagonal(H))⁻¹ Jᵀ W γ
update β .+= δ .* scales[perm]
"""


function mat_crunch(ψ,mat,q,wvs,puretor::Bool)::SparseMatrixCSC{Float64,Int}
#this function doesn't work but would be very helpful
   out = sparse(I,ψ.R.lng,ψ.R.lng)
   if puretor
      mat = sand(mat, wvs.ttp.vecs[:,:, ψ.σ] )
      mat = torsetter() # <------------
   end
   out = kron(tpart,out)
   elseif isnothing(wvs.ttp)
      out = kron(I(ψ.T.l), out)
   else # !isnothing(wvs.ttp.vals)
      out = kron(I(size(wvs.ttp.vals,1)),out)
   end #tor if
   #if
   #end vib
   droptol!(out,1e-11)
end

function sumder(rpid::Int,prm::Vector{Float64},ℋ::Vector{Op},ψ::Psi,wvs::Eigs,
                UR::SparseMatrixCSC{Float64,Int})::SparseMatrixCSC{Float64,Int}
   ind = rpid+1
   if ind ≤ length(ℋ) + HCLENGTH
      check = ℋ[ind-HCLENGTH].stg
      while check < zero(check)
         pm = prm[ind]
         out .+= enact(ℋ[ind],ψ,wvs, prm[ind-HCLENGTH, UR, true)
         ind += 1
         if ind-HCLENGTH ≤ length(ℋ)
            check = ℋ[ind-HCLENGTH].stg
         else
            check = 0
         end
      end
   end
   return out
end

function derivmat(rpid::Int,prm::Vector{Float64},ℋ::Vector{Op},ψ::Psi,wvs::Eigs,
                  UR::SparseMatrixCSC{Float64,Int})
   if ℋ[rpid].scl < 0 #should this be ≤ 0 ???
   elseif rpid ≤ 4 #pure rot
      temp = zeros(4)
      temp[rpid] = 1.0
      out = hrot2_hc(temp,ψ.R.N)
      out = kron() # <------------
   elseif 5 ≤ rpid ≤ 8 #spin-rot
      temp = zeros(4)
      temp[rpid-4] = 1.0
      out = hsr(temp,ψ.R)
      out = kron()  # <------------
   elseif 9 ≤ rpid ≤ 11 #qua
      temp = zeros(3)
      temp[rpid-8] = 1.0
      out = hqua(temp,ψ.R)
      out = kron()  # <------------
   elseif rpid ∈ hcount + 1 .+ 4*(0:XXX) # F
      q = floor(Int, 0.25*(rpid - 1))
      out = enact( Op(tf=OpFunc(Pt,2,q)), ψ,wvs,1.0, UR)
   elseif rpid ∈ hcount + 2 .+ 4*(0:XXX) # ρz
      q = floor(Int, 0.25*(rpid - 2))
      out = htsr2_hc(pr,wvs,ψ,σid)
   elseif rpid ∈ hcount + 3 .+ 4*(0:XXX) # ρx
      q = floor(Int, 0.25*(rpid - 3))
      out = htsr2_hc(pr,wvs,ψ,σid)
   elseif rpid ∈ hcount + 4 .+ 4*(0:XXX) # Vn
      q = floor(Int, 0.25*(rpid - 4))
      htor2_hc([zeros(3); 1.0], ψ.T.tops[q])
   else #user def
      out = enact(O::Op,ψ::Psi,wvs::Eigs, 1.0, UR, true)
      out .= sumder()
   end
   return tplus!(out)
end

function anaderiv(...)
   mat = derivmat(...)
   out = sand(mat, wvs.rst.vecs[] )
   return out #diag(out)
end

function jacob_term(...) # <----------
   ders = zeros(ψ.l,ψ.l,length(perm))
   for i ∈ 1:length(perm)
      ders[:,:,i] = derivmat(...) # <----------
   end
   for i ∈ 1:length(perm), j ∈ 1:length(perm)
      i==j ? J[jinds, i] : nothing
      ddmat[i,j, jinds] = der2_birss_elem(dx,dy,wvs,i) # <----------
   end
end

function der2_magnus_elem(dx,dy,h,e,v)
   2*sand(dx*pinv(e*I(size(h,1)) - h)*dy, v)
end
function der2_birss_elem(dx,dy,wvs,i)
   e_i = wvs.rst.vals[lims, i] # <----------
   vec_i = wvs.rst.vecs[lims, i] # <----------
   inds = filter(x->!isequal(x,i), lims)
   out = 0.0
   for j ∈ inds
      vec_j = wvs.rst.vecs[lims, j] # <----------
      e_j = wvs.rst.vals[lims, j] # <----------
      out += dx[i,j]*dy[j,i] / (e_i - e_j)
   end
   return 2.0*out
end



"""
Notes on an alternate cost function for robust fitting:
pseudo-Huber loss function:
L(x) = δ²( √(1 + ((k - f(x,y)+g(x,y))/δ)²) - 1)

d/dx L(x) = (-df/dx + dg/dx)*(k - f + g)/ ( √(1 + ((k - f(x,y)+g(x,y))/δ)²) )
"""
