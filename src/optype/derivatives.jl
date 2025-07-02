"""
This is the file that contains the analytic derivative functions
1st deriv - 1 stage - done!
1st deriv - 2 stage
jacobian - 1 stage - done!
jacobian - 2 stage
hessian - jtwj - done!
2nd deriv - 1 stage - done!
2nd deriv - 2 stage
hessian - birss 1 stage
hessian - birss 2 stage

20 may 25: code should be roughly finished. will fix from opt testing
"""

function sumder(out::SparseMatrixCSC{Float64,Int},
                id::Int,prm::Vector{Float64},stg::Vector{Int},
                ℋ::Vector{Op},ψ::Psi)
   ind = id+1
   if ind ≤ length(stg)+11
      check = stg[ind-11]
      while check < zero(check)
         pm = prm[ind]
         out .+= enact(ℋ[ind-11],ψ,pm)#,ur,ut)
         ind += 1
         if ind-11 ≤ length(stg)
            check = stg[ind-11]
         else
            check = 0
         end
      end
   end
   return out
end
function derivmat(id::Int,prm::Vector{Float64},scl::Vector{Float64},stg::Vector{Int},
                  ℋ::Vector{Op},ψ::Psi)
   if scl[id] ≤ 0 #should this be ≤ 0 ???
   elseif id ≤ 4 #pure rot
      pr = zeros(4)
      pr[id] = 1.0
      out = hrot2(pr,ψ.R)
      out = kron(I(ψ.T.lng),out)
   elseif 5 ≤ id ≤ 8 #spin-rot
      pr = zeros(5)
      pr[id-4] = 1.0
      out = hsr(pr,ψ.R)
      out = kron(I(ψ.T.lng),out)
   elseif 9 ≤ id ≤ 11 #qua
      pr = zeros(3)
      pr[id-9] = 1.0
      out = hqu(pr,ψ.R)
      out = kron(I(ψ.T.lng),out)
   else #user def
      #UR = ur(ψ.R.J,ψ.R.S)
      #UT = ur(Int( 0.5*(length(ψ.T.ms[1])-1) ))
      out = enact(ℋ[id-11],ψ,1.0)
      out .= sumder(out,id,prm,stg,ℋ,ψ)#,UR,UT)
   end
   return out
end
function anaderiv(id::Int,mc,prm::Vector{Float64},scl::Vector{Float64},stg::Vector{Int},
                  ℋ::Vector{Op},ψ::Psi,vec)

   mat = derivmat(id,prm,scl,stg,ℋ,ψ)
   U = sparse(ones(1))
   for i in 1:length(ψ.T.nf)
#having U at the end of kron maintians the tight block structure of the lower index tops
      if ψ.T.σ[i] == 0
         U = kron(ur(mc),U)
      else
         U = kron(sparse(1.0I,2mc+1,2mc+1),U)
   end;end
   U = kron(U,ur(ψ.R.J,ψ.R.S))
   mat = droptol!(U*mat*U,2*eps())
   mat = transpose(vec)*mat*vec
   return mat #diag(out)
end
function derivcalc(jlist,ℋ,ctrl,perm,vecs,prm,scl,stg)#::Matrix{Float64}
   s = ctrl.S
   mcalc = ctrl.mcalc
   msd = (2*mcalc+1)*length(ctrl.NFOLD)
   σs = σgen_indef(ctrl.NFOLD)
   σcnt = maximum(size(σs))
   derivs = zeros(Float64,size(vecs,2),σcnt,length(perm))
   for sc in 1:σcnt
      σ = σs[:,sc]
      ϕ = TPsi(ctrl.NFOLD,σ,mcalc)
      jsublist = jlist[isequal.(jlist[:,2],sc), 1] .* 0.5 #<-----NEED NEW JLIST FUNCTION
      for j in jsublist
         ψ = Psi(RPsi(j,s),ϕ)
         jd = Int(2.0*j) + 1
         dest = jvdest2(j,s,ctrl.vtmax) 
         #qns = qngen(j,s)
         vec = vecs[1:jd*msd,dest,sc]
         for i in 1:length(perm)
            #pid = perm[i]
            #ders = anaderive(pid,prm,scl,stg,ℋ,ψ,vec)
            #derivs[sind:find,sc,i] = ders#*scl[pid]
      derivs[dest,sc,i] = diag(anaderiv(perm[i],mcalc,prm,scl,stg,ℋ,ψ,vec))
         end#perm loop
      end #j loop
   end#σ loop
   return derivs
end#function
function build_jcbn!(jcbn,inds,jlist,ℋ,ctrl,perm,vecs,prm,scl,stg)
#   jcbn = zeros(Float64,size(inds,1),length(perm))
   deriv = derivcalc(jlist,ℋ,ctrl,perm,vecs,prm,scl,stg)
   for p in 1:length(perm)
      for a in 1:size(inds,1)
         #dν/dX = d/dX (ν_o - E_u + E_l) = -dE_u/dX + dE_l/dX 
         jcbn[a,p] = -deriv[inds[a,3],inds[a,2],p] + deriv[inds[a,6],inds[a,5],p]
      end
   end
   return jcbn
end

function sumder_2stg(out,j,s,nf,rpid,prm,stg,ops,ms,qns,tvecs)
   ind = rpid+1
   if ind ≤ length(stg)+18
      check = stg[ind-18]
      while check < zero(check)
         pm = prm[ind]
         out .+= enact_stg2(O,ψ,val,tvcs,mc,ur,ut)
         ind += 1
         if ind-18 ≤ length(stg)
            check = stg[ind-18]
         else
            check = 0
         end
      end
   end
   return out
end
function derivmat_2stg(j,s,nf,rpid,prm,scl,stg,ops,ms,qns,tvecs,mmax)
   if scl[id] ≤ 0 #should this be ≤ 0 ???
   elseif id ≤ 4 #pure rot
      pr = zeros(4)
      pr[id] = 1.0
      out = hrot2(pr,qns)
      out = kron(I(length(ms)),out)
   elseif 5 ≤ id ≤ 8 #spin-rot
      pr = zeros(5)
      pr[id-4] = 1.0
      out = hsr(pr,j,s,qns)
      out = kron(I(length(ms)),out)
   elseif 9 ≤ id ≤ 11 #qua
      pr = zeros(3)
      pr[id-9] = 1.0
      out = hqu(pr,j,s,qns)
      out = kron(I(length(ms)),out)
   else #user def
      out = enact_stg2(O,ψ,val,tvcs,mc,ur,ut)
      out .= sumder_2stg(out,id,prm,stg,ℋ,ψ)
   end
   return out
end
function anaderiv_2stg(prm,scl,stg,rpid,ops,j,s,nf,ms,qns,vec,tvec,mmax)
   mat = derivmat_2stg(j,s,nf,rpid,prm,scl,stg,ops,ms,qns,tvec,mmax)
   out = transpose(vec)*mat*vec
   return diag(out)
end
function derivcalc_2stg(jlist,ops,ctrl,perm,vecs,prm,scl,stg,tvecs)
   s = ctrl.S
   nf = ctrl.NFOLD
   σcnt = σcount(nf)
   derivs = zeros(Float64,size(vecs,2),σcnt,length(perm))
   msd = Int((ctrl.mmax+1)*(2s+1))
   for sc in 1:σcnt
      σ = sc - 1
      ms = msgen(nf,ctrl.mcalc,σ)
      tvcs = tvecs[:,:,sc]
      jsublist = jlist[isequal.(jlist[:,2],σ), 1] .* 0.5
      for j in jsublist
         #println(j)
         jd = Int(2.0*j) + 1
         sind, find = jvdest(j,s,ctrl.vtmax) 
         qns = qngen(j,s)
         vec = vecs[1:jd*msd,sind:find,sc]
         for i in 1:length(perm)
            pid = perm[i]
            ders = anaderiv_2stg(prm,scl,stg,pid,ops,j,s,nf,ms,qns,vec,tvcs,ctrl.mmax)
            derivs[sind:find,sc,i] = ders#*scl[pid]
         end#perm loop
      end #j loop
   end#σ loop
   return derivs
end#function



function build_hess!(hssn,jtw,jcbn,weights)
   mul!(jtw,jcbn',weights)
   mul!(hssn,jtw,jcbn)
   #return hssn, jtw
end
function build_hess(jtw,jcbn,weights)
   jtw = transpose(jcbn)*weights
   hssn = jtw*jcbn
   return hssn, jtw
end

function d2Ei_dxdy(vals,vecs,opx,opy,i,j,ψ)
   @assert i≠j
   veci = view(vecs,:,i)
   vecj = view(vecs,:,j)
   out = veci' * enact(opx,ψ) * vecj
   out *= vecj' * enact(opy,ψ) * veci
   out /= vals[i] - vals[j]
end

function invdiffmat(vals::Vector{Float64})::SparseMatrixCSC{Float64,Int}
   Δ = val .- val' # check prime ordering
   Δ = 1 ./Δ
   Δ[diagind(Δ)] .= 0
   return sparse(Δ)
end

"""
jh_levels calculates the jacobian & hessian for the /energy levels/
"""
function jh_levels(vals::Array{Float64},vecs::Array{Float64},prm::Vector{Float64},ℋ::Vector{Op})
   lp = length(perm)
   jac = zeros(length(vals,lp))
   hes = zeros(length(vals), Int( 0.5*lp*(1+lp) ))
   for j ∈ jσlist
      n = #size of basis set
      der_block = zeros(n,n,lp)
      for sc ∈ σcnt
         inds = dest = jvdest2(j,ctrl.S,ctrl.vtmax)
         vpart = vecs[1:n,inds,sc]
         Δi = invdiffmat(vals[inds,sc])
         for i ∈ lp
            der_block[:,:,i] = anaderiv()
            jac[inds,:] = diag(der_block[:,:,i])
         #end
         for j ∈ 1:i#i:length(perm)
            ind_der = (i-1)*lp + j
   hes[inds,inds,ind_der] = sum(der_block[:,:,i] .* Δ .* der_block[:,:,j],dims=2) 
         end;end#nested
      end#σs
   end#j
   return jac,hes
end#function 
function jh_levels_2stg(vals::Array{Float64},vecs::Array{Float64,2},
         tvecs::Array{Float64,2},prm::Vector{Float64},ℋ::Vector{Op})
   lp = length(perm)
   jac = zeros(length(vals,lp))
   hes = zeros(length(vals), Int( 0.5*lp*(1+lp) ))
   for j ∈ jσlist
      n = #size of basis set
      der_block = zeros(n,n,lp)
      for sc ∈ σcnt
         inds = dest = jvdest2(j,ctrl.S,ctrl.vtmax)
         vpart = vecs[1:n,inds,sc]
         Δi = invdiffmat(vals[inds,sc])
         for i ∈ lp
            der_block[:,:,i] = anaderiv_2stg()
            jac[inds,:] = diag(der_block[:,:,i])
         #end
         for j ∈ 1:i#i:length(perm)
            ind_der = (i-1)*lp + j
   hes[inds,inds,ind_der] = sum(der_block[:,:,i] .* Δ .* der_block[:,:,j],dims=2) 
         end;end#nested
      end#σs
   end#j
   return jac,hes
end#function 
function jh_build!(J,H,W,omc,vals,vecs,prm,ℋ,ctrl)
   jac,hes = jh_levels(vals,vecs,prm,ℋ)
   H .= 0.0
   for i ∈ eachindex(lins)
      u = inds[i,2]#???
      l = inds[i,4]#???
      J[i,:] .= jac[u,:] - jac[l,:]
      omc[i] .= lins[i] - vals[u] + vals[i]
      H .+= hess_reshape(hes[u,:] .- hes[l,:])*omc[i]*W[i]
   end
   H .+= J' * W * J
end
function jh_build!(J,H,W,omc,vals,vecs,prm,ℋ,ctrl)
   jac,hes = jh_levels_2stg(vals,vecs,tvecs,prm,ℋ)
   H .= 0.0
   for i ∈ eachindex(lins)
      u = inds[i,2]#???
      l = inds[i,4]#???
      J[i,:] .= jac[u,:] - jac[l,:]
      omc[i] .= lins[i] - vals[u] + vals[i]
      H .+= hess_reshape(hes[u,:] .- hes[l,:])*omc[i]*W[i]
   end
   H .+= J' * W * J
end

#these were stolen from the deprecated TriangularReshapes.jl
#and then some simple wrappers for shorter names
function vector_to_lower_triang!(M::AbstractMatrix{T}, v::AbstractVector{T}) where {T}
    n = size(M, 1)
    @assert size(M, 2) == n
    @assert length(v) >= n * (n + 1) / 2
    st = 1
    l = n
    Base.require_one_based_indexing(M, v)
    @inbounds @simd for i = 1:n
        for j = 0:(n - i)
            M[j + i, i] = v[st + j]
        end
        st += l
        l -= 1
    end
    return nothing
end
function lower_triang_to_vector!(v::AbstractVector{T}, M::AbstractMatrix{T}) where {T}
    n = size(M, 1)
    @assert size(M, 2) == n
    @assert length(v) >= n * (n + 1) / 2
    st = 1
    l = n
    Base.require_one_based_indexing(M, v)
    @inbounds @simd for i = 1:n
        for j = 0:(n - i)
            v[st + j] = M[j + i, i]
        end
        st += l
        l -= 1
    end
    return nothing
end
function v2lt!(M::AbstractMatrix{T}, v::AbstractVector{T}) where {T}
   vector_to_lower_triang!(M,v)
end
function v2lt(v::AbstractVector{T}) where {T}
   nM = Int(-0.5 + sqrt(0.25 + 2*length(v)))
   M = zeros(T, nM, nM)
   vector_to_lower_triang!(M,v)
   return M
end
function hess_reshape(v::AbstractVector{T}) where {T}
   nM = Int(-0.5 + sqrt(0.25 + 2*length(v)))
   M = zeros(T, nM, nM)
   vector_to_lower_triang!(M,v)
   return M - M'
end
function lt2v!(v::AbstractVector{T},M::AbstractMatrix{T}) where {T}
   lower_triang_to_vector!(v,M)
end
function lt2v!(v::AbstractVector{T},M::AbstractMatrix{T}) where {T}
   n = size(M,1)
   v = Vector{T}(undef, n * (n + 1) ÷ 2)
   lower_triang_to_vector!(v,M)
   return v
end

#=
given Bk = 1.75, Bn = 1.25, B± = 0.125, Dab = 0.2, N=1,S=0
prm = [1.75; 1.25; 0.125; 0.1]
val = [ 2.473791265187003
        4.026208734813001
        4.5]
v =[ -0.0918764   0.701112   0.707107;
     -0.991523   -0.129933   0.000000;
      0.0918764  -0.701112   0.707107]
zz = diagm([1;0;1])
tt = diagm(fill(2,3))
pm = rotl90(diagm([2;0;2]))
xz = √0.5 .* [0 -1 0; -1 0 1; 0 1 0]
Δ = val .- val'; Δ = 1 ./Δ; Δ[diagind(Δ)] .= 0; Δ
δ1 = [1/(val[1]-val[2]) 0; 0 1/(val[1]-val[3])]
δ2 = [1/(val[2]-val[1]) 0; 0 1/(val[2]-val[3])]
δ3 = [1/(val[3]-val[1]) 0; 0 1/(val[3]-val[2])]

=#












#yikes