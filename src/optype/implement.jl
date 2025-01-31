
"""
This file is the actual implementation of the information in types.jl & baseops.jl
"""

using DelimitedFiles

include("@__DIR__/../type.jl")
include("@__DIR__/../baseops.jl")
include("@__DIR__/../file_in.jl")

const csl = 29979.2458

function westereng(molnam::String)
   ctrl = ctrlinp(molnam)
   val, errs, ℋ, stgs = secordinp(molnam,ctrl)
   ℋ, stgs, errs = opreader(molnam,ctrl,ℋ,stgs,errs)
   ℋ = stgvalset(ℋ,stgs)
   #initialize vals
   #initialize vecs
   #initialize tvecs
   #initialize qns
for j in jmin:jmax
   Hrot = hrot2v2()
   if S≥1.0
      Hrot += Hsr + Hqua
   elseif S==0.5
      Hrot += Hsr
   end
   for σ in 1:σcnt
      ψ = Psi()
      vals[dest],vecs[dest] = tsrdiag_1(Hr,ℋ,ψ)
   end
end
end

#tsr_diag variants: 0 generic universal, 1 preloads Hrs, 2 two-stage preloads Hrs
function tsrdiag_0(ℋ::Vector{Op},ψ::Psi)
   H = enact(ℋ,ψ)
   U = sparse(ones(1))
   for i in 1:length(ψ.nf)
   if ψ.σ[i] == 0
      U = kron(u,ur(mcalc))
   else
      U = kron(u,sparse(1.0I,2mcalc+1,2mcalc+1))
   end;end
   U = kron(u,ur(ψ.J,ψ.S))
   H = U*H*U
   vals,vecs = eigen!(Symmetric(Matrix(H)))
   if (ctrl["assign"]=="ram36")||(ctrl["assign"]=="RAM36")
      perm = ramassign(vecs,j,s,mcalc,vtm)
      vals = vals[perm]
      vecs = vecs[:,perm]
   elseif ctrl["assign"]=="expectk"
      vals, vecs = expectkassign!(vals,vecs,j,s,nf,mcalc,σ)      
   elseif ctrl["assign"]=="eeo"
      vals, vecs = eeoassign!(vals,vecs,j,s,nf,mcalc,σ)
   else
      vals, vecs = expectassign!(vals,vecs,j,s,nf,mcalc,σ)
   end
   return vals, vecs
end

function tsrdiag_1(Hr::SparseMatrixCSC{Float64,Int},ℋ::Vector{Op},ψ::Psi)
   H = Hr + enact(ℋ,ψ)
   U = sparse(ones(1))
   for i in 1:length(ψ.nf)
   if ψ.σ[i] == 0
      U = kron(u,ur(mcalc))
   else
      U = kron(u,sparse(1.0I,2mcalc+1,2mcalc+1))
   end;end
   U = kron(u,ur(ψ.J,ψ.S))
   H = U*H*U
   vals,vecs = eigen!(Symmetric(Matrix(H)))
   if (ctrl["assign"]=="ram36")||(ctrl["assign"]=="RAM36")
      perm = ramassign(vecs,j,s,mcalc,vtm)
      vals = vals[perm]
      vecs = vecs[:,perm]
   elseif ctrl["assign"]=="expectk"
      vals, vecs = expectkassign!(vals,vecs,j,s,nf,mcalc,σ)      
   elseif ctrl["assign"]=="eeo"
      vals, vecs = eeoassign!(vals,vecs,j,s,nf,mcalc,σ)
   else
      vals, vecs = expectassign!(vals,vecs,j,s,nf,mcalc,σ)
   end
   return vals, vecs
end