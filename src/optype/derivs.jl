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


r_n = W/2 (k - ν_ij)² =  W/2 (k - E_i + E_j)²
χ² = γᵀ W γ

d r_n / dx = W(-dE_i/dx + dE_j/dx)*(k - E_i + E_j)
J_n = dE_i/dx - dE_j/dx
∇χ² = -Jᵀ W γ

d² r_n / dx² = W(-d²E_i/dx² + d²E_j/dx² )(k - E_i + E_j) + W(-dE_i/dx + dE_j/dx)²
d² r_n / dxdy = (-d²E_i/dxdy + d²E_j/dxdy )W(k - E_i + E_j) + W(-dE_i/dx + dE_j/dx)(-dE_i/dy + dE_j/dy)
S_xyn = d²E_i/dxdy - d²E_j/dxdy 
   = (2 ∑_k≠i ⟨i| dℋ/dx |k⟩⟨k| dℋ/dy |i⟩ / E_i - E_k') - (2 ∑_k≠j ⟨j| dℋ/dx |k⟩⟨k| dℋ/dy |j⟩ / E_j - E_k'))
H_xy = ∑_n d² r_n / dxdy 
H(χ²) = -∑_n S_xyn W_n γ_n + Jx_n W_n Jy_n

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
   out = sparse(I,ψ.R.lng,ψ.R.lng)
   if puretor
      if isone(length(ψ.nfs)) # cases 0,1,2
         #nothing
      elseif length(ψ.nfs)>1 && isnothing(wvs.top.vecs) # cases 3,4,5
         mat = torsetter!(ψ,O.q,mat)
      elseif length(ψ.nfs)>1 && !isnothing(wvs.top.vecs) && iszero(wvs.top.vecs[O.q].vecs)
         #nothing
      elseif length(ψ.nfs)>1 && !isnothing(wvs.top.vecs) && !iszero(wvs.top.vecs[O.q].vecs)
         mat = sand(mat,tvs[O.q].vecs[:,:, σ2ind(ψ.σs[O.q], O.q, ψ.nfs[O.q])]) #
         mat = torsetter!(ψ,O.q,mat)
      else
         @warn "unexpected condition in evalulation of tor op"
      end
   end # puretor
   out = kron(mat,out)
   elseif isnothing(wvs.ttp)
      out = kron(I(ψ.T.l), out)
   else # !isnothing(wvs.ttp.vals)
      out = kron(I(size(wvs.ttp.vals,1)),out)
   end #tor if
   return droptol!(out,1e-11)
end

function sumder(rpid::Int,prm::Vector{Float64},ℋ::Vector{Op},ψ::Psi,wvs::Eigs,
                UR::SparseMatrixCSC{Float64,Int})::SparseMatrixCSC{Float64,Int}
   ind = rpid+1
   if ind ≤ length(ℋ) + HCLENGTH
      check = ℋ[ind-HCLENGTH].stg
      while check < zero(check)
         pm = prm[ind]
         out .+= enact(ℋ[ind],ψ,wvs, prm[ind-HCLENGTH], UR, true)*prm[rpid]
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

function derivop(rpid::Int,prm::Vector{Float64},ℋ::Vector{Op},ψ::Psi,wvs::Eigs,
                  UR::SparseMatrixCSC{Float64,Int})
   if rpid > HCLENGTH
      if ℋ[rpid].scl < 0 #should this be ≤ 0 ???
      end
   elseif rpid ≤ 4 #pure rot
      temp = zeros(4)
      temp[rpid] = 1.0
      out = hrot2_hc(temp,ψ.R.N)
      out = kron(I(size(wvs.ttp.vals,1)),out)
   elseif 5 ≤ rpid ≤ 8 #spin-rot
      temp = zeros(4)
      temp[rpid-4] = 1.0
      out = hsr(temp,ψ.R)
      out = kron(I(size(wvs.ttp.vals,1)),out)
   elseif 9 ≤ rpid ≤ 11 #qua
      temp = zeros(3)
      temp[rpid-8] = 1.0
      out = hqua(temp,ψ.R)
      out = kron(I(size(wvs.ttp.vals,1)),out)
   elseif rpid ∈ hcount + 1 .+ 4*(0:XXX) # F
      q = floor(Int, 0.25*(rpid - 1))
      out = htor2_hc( [1.0; zeros(3)], ψ.T.tops[q] )
      out = mat_crunch(ψ,mat,q,wvs,true)
   elseif rpid ∈ hcount + 2 .+ 4*(0:XXX) # ρz
      q = floor(Int, 0.25*(rpid - 2))
      out = htsr2_hc(pr,wvs,ψ,σid)
      out = mat_crunch(ψ,mat,q,wvs,false)
   elseif rpid ∈ hcount + 3 .+ 4*(0:XXX) # ρx
      q = floor(Int, 0.25*(rpid - 3))
      out = htsr2_hc(pr,wvs,ψ,σid)
      out = mat_crunch(ψ,mat,q,wvs,false)
   elseif rpid ∈ hcount + 4 .+ 4*(0:XXX) # Vn
      q = floor(Int, 0.25*(rpid - 4))
      out = htor2_hc([zeros(3); 1.0], ψ.T.tops[q])
      out = mat_crunch(ψ,mat,q,wvs,true)
   else #user def
      out = enact(ℋ[rpid-HCLENGTH], ψ,wvs, 1.0, UR, true)
      out .= sumder(rpid,prm, ℋ,ψ,wvs, UR)
   end
   return tplus!(out)
end

function anaderiv(rpid::Int,prm::Vector{Float64},ℋ::Vector{Op},ψ::Psi,wvs::Eigs,
                  UR::SparseMatrixCSC{Float64,Int}, jinds::UnitRange{Int})
   mat = derivop(rpid,prm, ℋ,ψ,wvs, UR)
   out = sand(mat, wvs.rst.vecs[1:ψ.l,jinds,ψ.σ] )
   return droptol(sparse(out), 1e-10)
end

function jacob_term(perm, prm, ℋ,ψ,wvs, UR) # <----------
   jinds = jinds(ψ.R.J, ψ.R.S, ctrl.vtmax+1) 
   ders = zeros(ψ.l,ψ.l,length(perm))
   for i ∈ 1:length(perm)
      ders[:,:,i] = anaderiv(perm[i]Sp, prm, ℋ,ψ,wvs, UR, jinds)
   end
   return ders
end

function dEcalc(ctrl,prm,ℋ,wvs, perm, jσlist)
   σs = σgen(ctrl.NFOLD)
   X = dgen(ctrl.jmax)*dgen(ctrl.S)*(ctrl.vtmax+1)
   #out = zeros(X, size(wvs.rst.vals,1) ,length(perm),size(σs,2))
   J_eng = zeros( size(wvs.rst.vals,1), length(perm) )
   H_eng = zeros( length(perm), length(perm), size(wvs.rst.vals,1) ) 
   for i ∈ 1:size(jσlist,1)
      j,σ = jσlist[i,:]
      ψ = Psi( RPsi(j,ctrl.S), TTPsi(ctrl.NFOLD,σs[:,σ],ctrl.mcalc), σ )
      UR = ???? # <------
      inds = ???? # <------
      temp = jacob_term(perm, prm, ℋ,ψ,wvs, UR)
      J_eng[inds, :] = diag(temp)
      H_eng[:,:,inds] = der2_block()
   end
   return J_eng, H_eng
end

function der2_birss_elem(dx,dy,wvs,i, inds)
   e_i = wvs.rst.vals[inds, i]
   vec_i = wvs.rst.vecs[inds, i]
   inds = filter(x->!isequal(x,i), inds)
   out = 0.0
   for j ∈ inds
      vec_j = wvs.rst.vecs[inds, j]
      e_j = wvs.rst.vals[inds, j]
      out += dx[i,j]*dy[j,i] / (e_i - e_j)
   end
   return 2.0*out
end
function der2_block(ders,wvs,inds)
   temp = zero(ders)
   for i ∈ 1:length(inds), x ∈ 1:size(ders,3), y ∈ 1:size(ders,3)
      temp[x,y,i] = der2_birss_elem(ders[:,:,x],ders[:,:,y], wvs, i, inds)
   end
   return temp
end

function dE2dfconv!(Jf,Hf, Je,He, W,γ, linds, perm)
   #W is just inverse freq unc
   # γ = W * (ofreq .- cfreq)
   #Jf = zeros(size(linds,1), length(perm) )
   #Hf = zeros(length(perm), length(perm) )
   Jf .= (Je[linds[:,1],:] .- Je[linds[:,2],:]) .* W  # <------------------
   S = He[:,:,linds[:,1]] .- He[:,:,linds[:,2]]
   Hf .= -sum(x->S[:,:,x] * W[x] * γ[x], eachindex(γ)) + Jf' * Jf
   return Jf, Hf
end

function jac_hess_calc()
   J_eng, H_eng = dEcalc()
   J_frq, H_frq = dE2dfconv()
end


"""
function der2_magnus_elem(dx,dy,h,e,v)
   2*sand(dx*pinv(e*I(size(h,1)) - h)*dy, v)
end

Notes on an alternate cost function for robust fitting:
pseudo-Huber loss function:
L(x) = δ²( √(1 + ((k - f(x,y)+g(x,y))/δ)²) - 1)

d/dx L(x) = (-df/dx + dg/dx)*(k - f + g)/ ( √(1 + ((k - f(x,y)+g(x,y))/δ)²) )
"""
