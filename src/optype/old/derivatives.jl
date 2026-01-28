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
"""

function derivmat(id::Int,mc,prm::Vector{Float64},scl::Vector{Float64},stg::Vector{Int},
                  ℋ::Vector{Op},ψ::Psi)
   topcount = length(ψ.T.nf)
   if scl[id] < 0 #should this be ≤ 0 ???
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
   elseif 12 ≤ id ≤ 11+6*topcount
      #get fucked past wes
      opid = div(id-12, topcount)
      tid = mod(id, topcount) + 1
      if opid==0# F
#         out = kron(p_top(ψ.T, 2,tid),I(ψ.R.lng))
         out = htorhc(nf,1.0,0.0,mc,ψ.T.ms[tid],ψ.T.σ[tid])
         out = kron(out, I(ψ.R.lng))
      elseif opid==1# ρz
         out = kron(p_tor(ψ.T, 1,tid), Diagonal(nz(ψ.R.K, 1)))
      elseif opid==2# ρx
         out = kron(p_tor(ψ.T, 1,tid), nx(ψ.R, 1))
      elseif opid==3# V
         out = htorhc(ψ.T.nf[tid],0.0,1.0,mc,ψ.T.ms[tid],ψ.T.σ[tid])
         out = kron(out, I(ψ.R.lng))
      elseif opid==4# ηz
         out = kron(p_tor(ψ.T, 1,tid), sz(ψ.R, 1))
      elseif opid==5# ηx
         out = kron(p_tor(ψ.T, 1,tid), sx(ψ.R, 1))
      else #FUCK 
         @warn "something is very wrong with the torsional derivatives"
      end
      torsetter!(ψ.T,tid,out)

   else #user def
      out = enact(ℋ[id-11-6*topcount],ψ,1.0)
      #sumder!(out,id,prm,stg,ℋ,ψ)
   end
   return out
end
function derivmat_2stg(id::Int,mc,prm::Vector{Float64},scl::Vector{Float64},stg::Vector{Int},
                  ℋ::Vector{Op},ψ::Psi,tvecs::Array{Float64,2})
   topcount = length(ψ.T.nf)
   if scl[id] < 0 #should this be ≤ 0 ???
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
   elseif 12 ≤ id ≤ 11+6*topcount
      #get fucked past wes
      opid = div(id-12, topcount)
      tid = mod(id, topcount) + 1
      if opid==0# F
#         out = kron(p_top(ψ.T, 2,tid),I(ψ.R.lng))
         out = sand(htorhc(nf,1.0,0.0,mc,ψ.T.ms[tid],ψ.T.σ[tid]), tvecs)
         torsetter!(ψ.T,tid,out)
         out = kron(out, I(ψ.R.lng))
      elseif opid==1# ρz
         out = kron(sand(p_tor(ψ.T, 1,tid), tvecs), Diagonal(nz(ψ.R.K, 1)))
         torsetter!(ψ.T,tid,out)
      elseif opid==2# ρx
         out = kron(sand(p_tor(ψ.T, 1,tid), tvecs), nx(ψ.R, 1))
         torsetter!(ψ.T,tid,out)
      elseif opid==3# V
         out = sand(htorhc(ψ.T.nf[tid],0.0,1.0,mc,ψ.T.ms[tid],ψ.T.σ[tid]), tvecs)
         out = kron(out, I(ψ.R.lng))
         torsetter!(ψ.T,tid,out)
      elseif opid==4# ηz
         out = kron(sand(p_tor(ψ.T, 1,tid), tvecs), sz(ψ.R, 1))
         torsetter!(ψ.T,tid,out)
      elseif opid==5# ηx
         out = kron(sand(p_tor(ψ.T, 1,tid), tvecs), sx(ψ.R, 1))
         torsetter!(ψ.T,tid,out)
      else #FUCK 
         @warn "something is very wrong with the torsional derivatives"
      end
   else #user def
      out = enact_stg2(ℋ[id-11-6*topcount],ψ,1.0,tvecs)
      sumder_2stg!(out,id,prm,stg,ℋ,ψ,tvecs)
   end
   return out
end

function sumder!(out::SparseMatrixCSC{Float64,Int},
                id::Int,prm::Vector{Float64},stg::Vector{Int},
                ℋ::Vector{Op},ψ::Psi)::SparseMatrixCSC{Float64,Int}
   ind = id+1
   #tpcnt = length(ψ.T.nf)
   shift = 11+6*length(ψ.T.nf)
   if ind ≤ length(stg)
      check = stg[ind-shift]
      while check < zero(check)
         pm = prm[ind]
         out .+= enact(ℋ[ind-shift],ψ,pm)
         ind += 1
         if ind-shift ≤ length(stg)
            check = stg[ind-shift]
         else
            check = 0
         end#if
      end#while
   end#if
   #eturn out
end
function sumder_2stg!(out::SparseMatrixCSC{Float64,Int},
                id::Int,prm::Vector{Float64},stg::Vector{Int},
                ℋ::Vector{Op},ψ::Psi,tvecs)::SparseMatrixCSC{Float64,Int}
   ind = id+1
   tpcnt = length(ψ.T.nf)
   if ind ≤ length(stg)
      check = stg[ind-11-6*tpcnt]
      while check < zero(check)
         pm = prm[ind]
         out .+= out = enact_stg2(ℋ[id-11-6*topcount],ψ,1.0,tvecs)
         ind += 1
         if ind-11 ≤ length(stg)
            check = stg[ind-11]
         else
            check = 0
         end#if
      end#while
   end#if
end


function anaderiv(id::Int,mc,prm::Vector{Float64},scl::Vector{Float64},stg::Vector{Int},
                  ℋ::Vector{Op},ψ::Psi,vec)
   mat = derivmat(id,mc,prm,scl,stg,ℋ,ψ)
   U = kron(sparse(1.0I,ψ.T.lng,ψ.T.lng), ur(ψ.R.J,ψ.R.S))
   mat = droptol!(sand(mat,U),2*eps())

   mat = transpose(vec)*mat*vec
   return mat #diag(out)
end
function anaderiv_2stg(id::Int,mc,prm::Vector{Float64},scl::Vector{Float64},stg::Vector{Int},
                  ℋ::Vector{Op},ψ::Psi,vec,tvecs)
   mat = derivmat(id,mc,prm,scl,stg,ℋ,ψ,tvecs)
   U = kron(sparse(1.0I,ψ.T.lng,ψ.T.lng), ur(ψ.R.J,ψ.R.S))
   mat = droptol!(sand(mat,U),2*eps())

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
function jh_levels(vals::Array{Float64},vecs::Array{Float64},
                   prm::Vector{Float64},ℋ::Vector{Op})
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