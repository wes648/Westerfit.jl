
"""
This file is the actual implementation of the information in types.Jl & baseops.Jl
"""

using DelimitedFiles
using WIGXJPFjl

include("@__DIR__/../type.jl")
include("@__DIR__/../baseops.jl")
include("@__DIR__/../file_in.jl")
include("@__DIR__/../hc_ham.jl")
include("@__DIR__/../assign.jl")
include("@__DIR__/../ntop.jl")

const csl = 29979.2458
"""
powneg1(x) takes a number and returns (-1)^x. I realize this is a stupid looking
function to have but it evalutes every so slightly faster than just (-1)^x
"""
powneg1(k::Real)::Int = isodd(k) ? -1 : 1
δi(x::Real,y::Real)::Int = x==y

function srprep(J,S)
   ns = Δlist(J,S)
   nd = 2 .* Int.(ns) .+ 1
   ni = ones(Int, length(ns),2)
   ni[1,2] = nd[1]
   for i in 2:length(ns)
      ni[i,1] = ni[i-1,2] + 1
      ni[i,2] = ni[i,1] + nd[i] - 1
   end
   jd = Int((2.0*S+1.0)*(2.0*J+1.0))
   return ns, nd, ni, jd
end
function nindsgen(ns::UnitRange{Int})::Vector{UnitRange{Int}}
   ni = Vector{UnitRange{Int}}(undef,length(ns))
   ni[1] = 1:2*ns[1]+1
   for i in 2:length(ns)
      ni[i] = (ni[i-1][end]+1):(ni[i-1][end]+2ns[i]+1)
   end
   return ni
end

"""
This builds the rotational Wang Transformation matrix for every n in Δlist(j,s).
   This will be kronecker producted with an identy matrix of size 2*m+1 for the
   torsional-rotational nature. A purely rotational form can be built using m=0
"""
function ur(j::Float64,s::Float64)::SparseMatrixCSC{Float64, Int64}
   ns, nd, ni, jsd = srprep(j,s)
   out = spzeros(jsd,jsd)
   for i in 1:length(ns)
      r = ni[i,1]:ni[i,2]
      out[r, r] = ur(ns[i])
   end
   return out
end

function westereng(molnam::String)
   ctrl = ctrlinp(molnam)
   prm, errs, ℋ, stgs = secordinp(molnam,ctrl)
   ℋ, stgs, errs = opreader(molnam,ctrl,ℋ,stgs,errs)
   #ℋ = stgvalset(ℋ,stgs)
   σs = σgen_indef(ctrl["NFOLD"])
   σcnt = maximum(size(σs))
   sd = Int(2*ctrl["S"]+1)
   jlist = collect(0.5*iseven(sd):ctrl["Jmax"])
   mcd = Int(2*ctrl["mcalc"]+1)

   #initialize vals
   vtd = ctrl["vtmax"] + 1
   jfd = sd*Int(sum(2.0 .* jlist .+ 1.0))
   vals = zeros(Float64,jfd*vtd,σcnt)
   #initialize vecs
   if ctrl["stages"]==1
      vecs = zeros(Float64,Int(sd*(2*ctrl["Jmax"]+1)*mcd),jfd*vtd,σcnt)
   elseif ctrl["stages"]==2
      vl = sd*(2*ctrl["Jmax"]+1)*(ctrl["mmax"]+1)
      vecs = zeros(Float64,Int(vl),jfd*vtd,σcnt)
      #initialize tvecs
      tvals = zeros(ctrl["mmax"]+1,σcnt)
      tvecs = zeros(2*ctrl["mcalc"]+1,ctrl["mmax"]+1,σcnt)
   end
   #do the energy level calculation
   if ctrl["stages"]==1
      vals,vecs = tsrcalc_1stg!(vals,vecs,jlist,σs,ctrl,prm,stgs,ℋ)
   elseif ctrl["stages"]==2
   vals,vecs,tvals,tvecs = tsrcalc_2stg!(vals,vecs,tvals,tvecs,jlist,σs,ctrl,prm,stgs,ℋ)
   else
      @warn "Invalid stages number"
   end
   #qns = bigqngenv2()
   return vals,vecs #,qns
end

function jvdest2(j,s,vtm)
"""
This returns a unit range spanning from the first to the final indices for a 
   certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return snd:fnd
end

#tsrcalc2(prm,stg,cdo,nf,ctrl,jlist)
function tsrcalc_1stg!(vals,vecs,jlist,σs,ctrl,prm,stg,ℋ)
   σcnt = size(σs,2)
   msd = (2*ctrl["mcalc"]+1)*length(ctrl["NFOLD"])
   for j in jlist
      jd = Int(2j+1)
      dest = jvdest2(j,ctrl["S"],ctrl["vtmax"]) 
      ψ = Psi(j,ctrl["S"])
      Hrot = hrot2(prm[1:4],ψ)
      if ctrl["S"]≥1.0
         Hrot += Hsr(prm[5:8],ψ.J,ψ.S,ψ) + Hqua(prm[9:11],ψ.J,ψ.S,ψ)
      elseif ctrl["S"]==0.5
         Hrot += Hsr(prm[5:8],ψ.J,ψ.S,ψ)
      end
      Hrot = sparse(Symmetric(Hrot,:L))
      for sc in 1:σcnt
         ψ = Psi(j,ctrl["S"],nf=ctrl["NFOLD"],σ=σs[sc],mc=ctrl["mcalc"])
         vals[dest,sc],vecs[1:jd*msd,dest,sc] = tsrdiag_1(Hrot,ctrl,ℋ,ψ,σs[:,sc])
      end#σs
   end#j
   return vals, vecs
end#f

function torcalc!(tvals,tvecs,ctrl,ℋ,ϕ,stg,sc)
   tsize = (2*ctrl["mcalc"]+1)^length(ctrl["NFOLD"])
   dest = 1:ctrl["mmax"]+1
   H = torbuild(ℋ,ϕ,stg,tsize,ctrl["mcalc"])
   H = Matrix(H)
   H = eigen!(H)
   tvals[dest,sc] = H.values[dest]
   tvecs[:,dest,sc] = H.vectors[:,dest]
end

function tsrcalc_2stg!(vals,vecs,tvals,tvecs,jlist,σs,ctrl,prm,stg,ℋ)
   σcnt = size(σs,2)#this doesn't work for 1 top
   tsize = (2*ctrl["mcalc"]+1)^length(ctrl["NFOLD"])
   for sc in 1:σcnt
      ϕ = Psi(ctrl["NFOLD"],σs[sc],ctrl["mcalc"])
      torcalc!(tvals,tvecs,ctrl,ℋ,ϕ,stg,sc)
      #tvals,tvecs = eigen!(Matrix(torbuild(ℋ,ϕ,stg,tsize)))
   end
   msd = Int(2ctrl["S"]+1)*(ctrl["mmax"]+1)
for j in jlist
   jd = Int(2j+1)
   dest = jvdest2(j,ctrl["S"],ctrl["vtmax"]) 
   ψ = Psi(j,ctrl["S"])
   Hrot = hrot2(prm[1:4],ψ)
   if ctrl["S"]≥1.0
      Hrot += Hqua(prm[9:11],ψ.J,ψ.S,ψ)
      if norm(prm[5:8]) > 0.0
         Hrot += Hsr(prm[5:8],ψ.J,ψ.S,ψ)
      end
   elseif ctrl["S"]==0.5
      Hrot += Hsr(prm[5:8],ψ.J,ψ.S,ψ)
   end
   Hrot = sparse(Symmetric(Hrot,:L))
   for sc in 1:σcnt
      ψ = Psi(j,ctrl["S"],nf=ctrl["NFOLD"],σ=σs[sc],mc=ctrl["mcalc"])
      vals[dest,sc],vecs[1:jd*msd,dest,sc] = tsrdiag_2(Hrot,ctrl,tvals[:,sc],tvecs[:,:,sc],ℋ,ψ,stg)
   end#σs
end#j
   return vals,vecs,tvals,tvecs
end#f

#tsr_diag variants: 0 generic universal, 1 preloads Hrs, 2 two-stage preloads Hrs
function tsrdiag_0(ℋ::Vector{Op},ψ::Psi)
   H = enact(ℋ,ψ)
   U = sparse(ones(1))
   for i in 1:length(ψ.nf)
#having U at the end of kron maintians the tight block structure of the lower index tops
   if ψ.σ[i] == 0
      U = kron(ur(mcalc),U)
   else
      U = kron(sparse(1.0I,2mcalc+1,2mcalc+1),U)
   end;end
   U = kron(u,ur(ψ.J,ψ.S))
   H = droptol!(U*H*U,2*eps())
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

function tsrdiag_1(Hr::SparseMatrixCSC{Float64,Int},ctrl,ℋ::Vector{Op},ψ::Psi,σs)
   H = kron(I(2*ctrl["mcalc"]+1)^length(ctrl["NFOLD"]),Hr) 
   #@show size(H)
   #@show size(enact(ℋ,ψ))
   H += enact(ℋ,ψ)
   U = sparse(ones(1))
   for i in 1:length(ψ.nf)
   if σs[i] == 0
      U = kron(U,ur(ctrl["mcalc"]))
   else
      U = kron(U,sparse(1.0I,2ctrl["mcalc"]+1,2ctrl["mcalc"]+1))
   end;end
   U = kron(U,ur(ψ.J,ψ.S))
   H = U*H*U
   vals,vecs = eigen!(Symmetric(Matrix(H)))
   if (ctrl["assign"]=="ram36")||(ctrl["assign"]=="RAM36")
      perm = ramassign(vecs,ψ.J,ψ.S,ctrl["mcalc"],ctrl["vtmax"])
      vals = vals[perm]
      vecs = vecs[:,perm]
   elseif ctrl["assign"]=="expectk"
      vals, vecs = expectkassign!(vals,vecs,ψ.J,ψ.S,nf,ctrl["mcalc"],σ)      
   elseif ctrl["assign"]=="eeo"
      vals, vecs = eeoassign!(vals,vecs,ψ.J,s,nf,mcalc,σ)
   else
      vals, vecs = expectassign!(vals,vecs,ψ.J,s,nf,mcalc,σ)
   end
   return vals, vecs
end
function tsrdiag_2(Hr::SparseMatrixCSC{Float64,Int},ctrl,tvals,tvecs,ℋ::Vector{Op},ψ::Psi,stg)
   #printstyled("ψ.J = $(ψ.J), ψ.σ = $(ψ.σ)\n",color=:cyan)
   H = kron(I(ctrl["mmax"]+1)^length(ctrl["NFOLD"]),Hr) 
   H[diagind(H)] += kron(tvals, ones(Int((2ψ.J+1)*(2ψ.S+1)) ))
   h_stg2build!(H,ℋ,ψ,stg,(2*ctrl["mcalc"]+1)^length(ctrl["NFOLD"]),tvecs,ctrl["mcalc"])
   U = sparse(ones(1))
   for i in 1:length(ψ.nf)
      l = ctrl["mmax"]+1
      U = kron(U,sparse(1.0I,l,l))
   end
   U = kron(U,ur(ψ.J,ψ.S))
   #@show H
   H = U*H*U
   vals,vecs = eigen!(Symmetric(Matrix(H),:L))
   #@show vecs
   perm = twostg_assign(vecs,ψ.J,ψ.S,ctrl["mmax"],ctrl["vtmax"])
   vals = vals[perm]
   vecs = vecs[:,perm]
   #@show nnz(dropzeros(sparse(vals)))
   return vals, vecs
end

function qnlabv2(j,s,vtm,σ)::Array{Int,2}
   nlist = Δlist(j,s)
   jsd = Int((2*j+1)*(2*s+1))
   vd = Int(vtm+1)
   out = zeros(Int,0,3)
   for n in nlist
      nd = Int(2*n+1)
      part = zeros(Int,nd,3)
      part[:,1] = fill(n,nd)
      part[:,2] = collect(Int,-n:n)
      part[:,3] = k2kc.(part[:,1],part[:,2])
      out = vcat(out,part)
   end
   out[:,2] = abs.(out[:,2])
   out = kron(ones(Int,vd),out)
   vtrray = kron(collect(0:vtm) ,ones(Int,jsd))
   out = hcat(fill(Int(2*j),size(out,1)),out,vtrray,fill(σ,jsd*vd))
   return out
end
