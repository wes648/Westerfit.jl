"""
This is the file that contains the analytic derivative functions
jacobian - 1 stage - done!
jacobian - 2 stage
hessian - jtwj - done!
hessian - briss 1 stage
hessian - briss 2 stage
"""

function sumder(out::SparseMatrixCSC{Float64,Int},
                id::Int,prm::Vector{Float64},stg::Vector{Int},
                ℋ::Vector{Op},ψ::Psi)
   ind = id+1
   if ind ≤ length(stg)+11
      check = stg[ind-11]
      while check < zero(check)
         pm = prm[ind]
         out .+= enact(ℋ[ind-11],ψ,pm,ur,ut)
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
      out = enact(ℋ[id-11],ψ,prm[id],ur,ut)
      out .= sumder(out,id,prm,stg,ℋ,ψ)
   end
   return out
end

function anaderiv(id::Int,prm::Vector{Float64},scl::Vector{Float64},stg::Vector{Int},
                  ℋ::Vector{Op},ψ::Psi,vec::Vector{Float64})::Vector{Float64}
   mat = derivmat(id,prm,scl,stg,ℋ,ψ)
   out = transpose(vec)*mat*vec
   return diag(out)
end

function derivcalc(jlist,ℋ,ctrl,perm,vecs,prm,scl,stg)::Matrix{Float64}
   s = ctrl["S"]
   mcalc = ctrl["mcalc"]
   σs = σgen_indef(ctrl["NFOLD"])
   σcnt = maximum(size(σs))
   derivs = zeros(Float64,size(vecs,2),σcnt,length(perm))
   for sc in 1:σcnt
      σ = σs[:,sc]
      ϕ = TPsi(ctrl["NFOLD"],σ,mcalc)
      jsublist = jlist[isequal.(jlist[:,2],sc-1), 1] .* 0.5 #<-----NEED NEW JLIST FUNCTION
      for j in jsublist
         ψ = Psi(RPsi(j,s),ϕ)
         jd = Int(2.0*j) + 1
         sind, find = jvdest(j,s,ctrl["vtmax"]) 
         qns = qngen(j,s)
         vec = vecs[1:jd*msd,sind:find,sc]
         for i in 1:length(perm)
            #pid = perm[i]
            #ders = anaderive(pid,prm,scl,stg,ℋ,ψ,vec)
            #derivs[sind:find,sc,i] = ders#*scl[pid]
      derivs[sind:find,sc,i] = anaderive(perm[i],prm,scl,stg,ℋ,ψ,vec)
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
