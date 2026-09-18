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

nfit_count(ℋ::Vector{Term})::Int = sum(map(x->!iszero(x.scl), ℋ))

function derivop_0(T::Term, ψ::Psi, wvs::Eigs,
                 UR::SparseMatrixCSC{Float64,Int})::SparseMatrixCSC{NUMTYPE,Int}
   out = enact(T.scl, T.ops[1], ψ,wvs, UR)
   @inbounds for i ∈ 2:T.l
      out += enact(T.scl, T.ops[i], ψ,wvs, UR)
   end
   return droptol!( tplus!(out), 1e-11)
end

function anaderiv(T::Term,ψ::Psi,wvs::Eigs,
                  UR::SparseMatrixCSC{Float64,Int}, jinds::UnitRange{Int})
   mat = derivop_0(T, ψ,wvs, UR)
   L = 1:size(mat,1)
   #@show ψ.R.J
   #@show ψ.σ
   out = sand(mat, wvs.rst.vecs[ L , jinds, ψ.σ] )
   return droptol!(sparse(out), 1e-10)
end

function jacob_term(ctrl,ℋ::Vector{Term},ψ::Psi,wvs::Eigs, UR)
   jind = jinds(ψ.R.J, ψ.R.S, ctrl.vtmax+1) 
   L = (ψ.R.lng * (ctrl.vtmax + 1))
   ders = zeros(L,L, nfit_count(ℋ))
   j = 1
   for i ∈ eachindex(ℋ)
      if !iszero(ℋ[i].scl)
         ders[:,:,j] = anaderiv(ℋ[i], ψ,wvs, UR, jind)
         j += 1
      end
   end
   return ders
end

function dEcalc(ctrl,ℋ,wvs, jσlist)
   σs = σgen(ctrl.NFOLD)
   J_eng = zeros( size(wvs.rst.vals,1), size(σs,2), nfit_count(ℋ) )
   nprm = nfit_count(ℋ)
   H_eng = zeros( nprm, nprm, size(wvs.rst.vals,1), size(σs,2) ) 
   #@show size(J_eng)
   for i ∈ 1:size(jσlist,1)
      j,σ = jσlist[i,:]
      j *= 0.5
   #   @show σ
      ψ = Psi( RPsi(j,ctrl.S), TTPsi(ctrl.NFOLD,σs[:,σ],ctrl.mcalc), σ )
      UR = ur(ψ.R.J, ψ.R.S)
      inds = jinds(ψ.R.J, ψ.R.S, ctrl.vtmax+1) 
      temp = jacob_term(ctrl,ℋ,ψ,wvs, UR)
      #@show size(temp)
      for l ∈ 1:nprm
         J_eng[inds, σ, l] = diag(temp[:,:,l])
      end
      if ctrl.trueHess
         H_eng[:,:,inds,σ ] = der2_block(temp,wvs,inds, σ)
      end
   end
   return J_eng, H_eng
end

function der2_birss_elem(dx,dy,wvs,i, inds, σ)
   e_i = wvs.rst.vals[inds[i], σ]
   #vec_i = wvs.rst.vecs[inds, i]
   shift = -minimum(inds) + 1
   tinds = filter(x->!isequal(x,inds[i]), inds)
   out = 0.0
   for j ∈ tinds
      #vec_j = wvs.rst.vecs[inds, j]
      out += dx[i, j+shift]*dy[j+shift, i] / (e_i - wvs.rst.vals[j, σ])
   end
   return 2.0*out
end
function der2_block(ders,wvs,inds,σ)
   temp = zeros(size(ders,3), size(ders,3), length(inds))
#   @show size(ders)
#   @show size(ders,3)
#   @show length(inds)
   for i ∈ 1:length(inds), x ∈ 1:size(ders,3), y ∈ 1:size(ders,3)
#      println("i = $i, x = $x, y = $y")
      temp[x,y,i] = der2_birss_elem(ders[:,:,x],ders[:,:,y], wvs, i, inds, σ)
   end
   return temp
end

function dE2dfconv!(flag::Bool,Jf,Hf, Je,He, W,γ, linds)
   #W is just inverse freq unc
   # γ = W * (ofreq .- cfreq)
   Hf .= zero(Hf)
   for σ ∈ unique(linds[:,3])
      dest = findall(x->x==σ, linds[:,3])
      uσl = linds[dest, 2]
      lσl = linds[findall(x->x==σ, linds[:,6]), 5]
      Jf[dest,:] .= (Je[uσl, σ+1, :] .- Je[lσl, σ+1, :]) .* W[dest]
      if flag
         S = He[:,:,uσl,σ+1] .- He[:,:,lσl,σ+1]
         Hf .+= sum(x->S[:,:,x] * W[dest[x]] * γ[dest[x]], eachindex(γ[dest]))
      end
   end
   Hf .+= Jf' * Jf
   #@show size(Hf)
   return Jf, Hf
end

function deriv_calc!(Jf,Hf, ctrl,ℋ,wvs,lins, γ)
   jσlst = jlister(lins.inds)
   Je, He = dEcalc(ctrl, ℋ,wvs, jσlst)
   dE2dfconv!(ctrl.trueHess, Jf,Hf, Je,He, lins.wght,γ, lins.inds)
   return Jf, Hf
end
function deriv_calc(ctrl,ℋ,wvs,lins, γ)
   nprm = nfit_count(ℋ)
   Jf = zeros(length(γ), nprm)
   Hf = zeros(nprm,nprm)
   jσlst = jlister(lins.inds)
   @time "dE/dx" Je, He = dEcalc(ctrl, ℋ,wvs, jσlst)
   @time "dν/dx" dE2dfconv!(ctrl.trueHess,Jf,Hf, Je,He, lins.wght,γ, lins.inds)
   if iszero(Jf)
      println("FUCK JACOBIAN IS ZERO")
   end
   if iszero(Hf)
      println("FUCK HESSIAN IS ZERO")
   end
   return Jf, Hf
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
