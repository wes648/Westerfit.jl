
#using LinearAlgebra, SparseArrays, WIGXJPFjl, BenchmarkTools
"""
eh(N) returns the square root of the N²|NK⟩ matrix element, √(N(N+1)). See eh2
for the non-square rooted version.
"""
eh(x::Real)::Float64 = √(x*(x+1))
"""
eh2(N) returns the the N²|NK⟩ matrix element, N(N+1). See eh for the
automatically square rooted version.
"""
eh2(x::Real)::Float64 = x*(x+1)
□rt(x::Real)::Float64 =√(x*(x>zero(x)))
fh(x::Real,y::Real)::Float64 = □rt((x-y)*(x+y+1))
jnred(j::Real,n::Real)::Float64 = √((2*j+1)*(2*n+1))
nred(n::Real)::Float64 = √(n*(n+1)*(2*n+1))
"""
powneg1(x) takes a number and returns (-1)^x. I realize this is a stupid looking
function to have but it evalutes every so slightly faster than just (-1)^x
"""
powneg1(k::Real)::Int = isodd(k) ? -1 : 1
"""
δ(x,y) takes two number and returns the Kronecker delta as a float. See δi for
the integer version
"""
δ(x::Real,y::Real)::Float64 = x==y
"""
δi(x,y) takes two number and returns the Kronecker delta as an integer. See δ
for the float version
"""
δi(x::Real,y::Real)::Int = x==y
T(l::Int,q::Int)::Int = l*(l+1) + q + 1
Tq(q::Int)::Int = 3 + q #quadrupole variant (only has 2nd rank components)
Tsr(l::Int,q::Int)::Int = δi(l,2) + abs(q) + 1 #Cs sr version, no 1st rk, & symm

"""
tplus!(a) replaces the matrix a with the sum of it and it's transpose. A dense
and sparse variant are available
"""
tplus!(a::Array{Float64,2})::Array{Float64,2} = a .+= permutedims(a)
function tplus!(a::SparseMatrixCSC{Float64, Int64})::SparseMatrixCSC{Float64, Int64}
   a .+= permutedims(a)
end

"""
qngen(j,s) generates the quantum numbers that the J dependent parts of the 
Hamiltonian processes. J is the total angular moment and S is the spin.
Returns a 2D array with Ns in the first column and Ks in the second 
"""
function qngen(j,s)::Array{Int,2}
   ns, nd, ni, jsd = srprep(j,s)
   out = zeros(Int,jsd,2)
   for i in 1:length(ns)
      out[ni[i,1]:ni[i,2],1] .= ns[i]
      out[ni[i,1]:ni[i,2],2] = collect(Int,-ns[i]:ns[i])
   end
   #[n k]
   return out
end
function qnlabv(j,s,nf,vtm,σ)::Array{Int,2}
   σt = σtype(nf,σ)
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
   vtrray = kron(nf .* vtcoll(vtm,σt) .+ σ,ones(Int,jsd))
   out = hcat(fill(Int(2*j),size(out,1)),out,vtrray,fill(σ,jsd*vd))
   return out
end
####   new 2nd order   ####
hr2on(ns,ks,bk,bn) = @. bn*eh2(ns) + bk*ks^2 
hr2of1(ns,ks,dab) = @. dab*(ks-0.5)*fh(ns,ks-1)
hr2of2(ns,ks,bpm) = @. bpm*fh(ns,ks-1)*fh(ns,ks-2)
function hrot2(pr,qns)::SparseMatrixCSC{Float64, Int64}
   ns = view(qns,:,1)
   ks = view(qns,:,2)
   out = spzeros(size(ns,1),size(ns,1))
   #p0 = hr2on(ns,ks,pr[1],pr[2])
   out[diagind(out)] .= hr2on(ns,ks,pr[1],pr[2])
   #p1 = hr2of1(ns[2:end],ks[2:end], pr[4])
   out[diagind(out,1)] .= hr2of1(ns[2:end],ks[2:end], pr[4])
   #p2 = hr2of2(ns[3:end],ks[3:end], pr[3])
   out[diagind(out,2)] .= hr2of2(ns[3:end],ks[3:end], pr[3])
   #out = spdiagm(0=>p0,1=>p1,2=>p2)
   return Symmetric(dropzeros!(out))
end
"""
hrotest(pr,j,s) generates the 2nd order rotational Hamiltonian for the given J
and S pair. pr is an array of length 4 with values of BK, BN, B±, and Dab 
respectively
"""
function hrotest(pr,j,s)
   qns = qngen(j,s)
   out = hrot2(pr,qns)
   return out
end
function hrotest(pr,n)
   qns = qngen(n,0)
   out = hrot2(pr,qns)
   return out
end

function nsred(l::Int,nb::Int,nk::Int)::Float64
   return 0.5*( 
   √(nk*(nk+1)*(2*nk+1))*
      wig6j( 1, 1, l,
            nb,nk,nk)*powneg1(l) + 
   √(nb*(nb+1)*(2*nb+1))*
      wig6j( 1, 1, l,
            nk,nb,nb))
end
function jsred(j,s,nb::Int,nk::Int)::Float64
   return wig6j(nk, s, j,
                 s,nb, 1)*√((2*nb+1)*(2*nk+1))
end
function srelem(pr::Float64,l::Int,q::Int,j,s,nb,kb,nk,kk)::Float64
   return pr*wig3j( nb,l,nk,
                   -kb,q,kk)*√(2*l+1)*
       nsred(l,nb,nk)*jsred(j,s,nb,nk)*powneg1(nb-nk-kb)
end
const ts = SA[0  1 1 2  2 2  2 2;
              0 -1 1 0 -1 1 -2 2;
              1  2 2 3  4 4  5 5;
              1  1 1 1 -1 1  1 1]
function hsr(pr::Array{Float64},j,s,qns::Array{Int})::SparseMatrixCSC{Float64, Int64}
   ns = view(qns,:,1)
   ks = view(qns,:,2)
   le = size(ns,1)
   out = spzeros(le,le)
   #awkward special array of rank, component, prm ind, & sign
   #each col is a different parameter
   #ts = SA[0  1 1 2  2 2  2 2;
   #        0 -1 1 0 -1 1 -2 2;
   #        1  2 2 3  4 4  5 5;
   #        1  1 1 1 -1 1  1 1]
   for i in 1:8
      tv = ts[:,i]
      prm = pr[tv[3]]*tv[4]
      if prm ≠ 0.0
      for a in 1:le, b in 1:le
         nb = ns[b]
         kb = ks[b]
         nk = ns[a]
         kk = ks[a]
         if abs(nb-nk)≤1 && (tv[2]+kk-kb)==0
            out[b,a] += srelem(prm,tv[1],tv[2], j,s,nb,kb,nk,kk)
         end#selection rule if
      end#sr ind for loop
      end#prm chck if
   end#sr term for loop
   dropzeros!(out)
   out .*= nred(s)*powneg1(j+s)
   return out
end

function wiginv(s::Real)::Float64
   if s<one(s)
      return 0.0
   else
      return inv(wig3j( s,2,s,
                       -s,0,s))
   end
end
function qured(j,s,nb,nk)
   return jnred(nb,nk)*
          wig6j(j, s,nb,
                2,nk, s)
end
function qulm(pr,q,j,s,nb,kb,nk,kk)#::Array{Float64,2}
   return pr*qured(j,s,nb,nk)*
#             δ(nb,nk)* #This line can be used to emulate the perturbative 
             wig3j( nb, 2,nk,
                   -kb, q,kk)*powneg1(nb+nk-kb)
end
#const tq = SA[0 -1 1 -2 2; 
#              1  2 2  3 3;
#              1 -1 1  1 1]
function hqu(pr,j,s,qns)::SparseMatrixCSC{Float64, Int64}
   ns = view(qns,:,1)
   ks = view(qns,:,2)
   le = size(ns,1)
   out = spzeros(le,le)
   #awkward special array of rank, component, prm ind, & sign
   #each col is a different parameter
   tq = SA[0 -1 1 -2 2; 
           1  2 2  3 3;
           1 -1 1  1 1]
   for i in 1:5
      tv = tq[:,i]
      prm = pr[tv[2]]*tv[3]
      if prm ≠ 0.0
      for a in 1:le, b in 1:le
         nb = ns[b]
         kb = ks[b]
         nk = ns[a]
         kk = ks[a]
         if abs(nb-nk)≤2 && (tv[1]+kk-kb)==0
            out[b,a] += qulm(prm,tv[1], j,s,nb,kb,nk,kk)
         end#selection rule if
      end#qu ind for loop
      end#prm chck if
   end#qu term for loop
   dropzeros!(out)
   out .*= 0.25*wiginv(s)*powneg1(j+s)
   #@show out
   return out
end

function htor2(sof::Array{Float64},ms::Array{Int})::SparseMatrixCSC{Float64, Int64}
   out = sof[1]*pa_op(ms,2)
   out += sof[4].*(I(size(out,1)) .- cos_op(ms,1))
   return out
end
function htor2v2(sof::Array{Float64},nf::Int,mc::Int,σ::Int)
   if nf≠0
      ms = msgen(nf,mc,σ)
      out = sof[1]*pa_op(ms,2)
      out += sof[4].*(I(size(out,1)) .- cos_op(ms,1))
   else
      out = [0.0]
      ms = [1]
   end
   return out, ms
end

####   individual operators   ####

function nnss_check(a,b)::Int
   a = a*iseven(a) + (a-1)*isodd(a)
   b = b*iseven(b) + (b-1)*isodd(b)
   return min(a,b)
end
ns_el(j,s,p,n)::Float64 = (0.5*eh2(j) - eh2(n) - eh2(s))^p
function nnss_op(j,s,qns,a,b)::Diagonal{Float64, Vector{Float64}}
   if (a≠0)&&(b≠0)
      c = nnss_check(a,b)
      a -= c
      b -= c
      @views out = eh.(qns[:,1]).^a .* ns_el.(j,s,c,qns[:,1]) .* eh(s)^b
   elseif (a≠0)&&(b==0)
      @views out = eh.(qns[:,1]).^a
   elseif (a==0)&&(b≠0)
      out = fill(eh(s)^b,size(qns,1))
   else
      out = ones(size(qns,1))
   end
   return Diagonal(out)
end
nt2_op(qns,p)::Diagonal{Float64, Vector{Float64}} = @views out = Diagonal(eh2.(qns[:,1]).^p)
function nz_op(qns,p)::Diagonal{Float64, Vector{Float64}} 
   if p≠0
      @views out = Diagonal(qns[:,2].^p)
   else
      out = Diagonal(ones(size(qns,1)))
   end
   return out
end

function np_op(qns,p::Int)::SparseMatrixCSC{Float64, Int64}
   ns = qns[1+p:end,1]
   part = ones(length(ns))
   if p ≤ length(ns) && p ≠ 0
      ks = qns[1+p:end,2]
      part = ones(length(ks))
      for o in 1:p
         part .*= fh.(ns,ks.-o)
      end
   #end#original if
      out = spzeros(size(qns,1),size(qns,1))
      out[diagind(out,-p)] = part
   elseif p > length(ns) && p ≠ 0
      out = spzeros(size(qns,1),size(qns,1))
      out[diagind(out,-p)] = part
   else
      out = spdiagm(ones(size(qns,1)))
   end
   return out
end
npm_op(qns::Matrix{Int64},p::Int) = Symmetric(np_op(qns,p),:L)
function iny_op(qns::Matrix{Int64},p::Int)::SparseMatrixCSC{Float64,Int}
   if p≠0
      out = np_op(qns,1-δi(p,0))
      out .-= permutedims(out)
      return dropzeros!(out)
   else
      return spdiagm(ones(size(qns,1)))
   end
end

function sqpart(j,s,q,bqn,kqn)::Float64
   nb = bqn[1]
   kb = bqn[2]
   nk = kqn[1]
   kk = kqn[2]
   return wig3j(nb,1,nk,-kb,q,kk)*wig6j(s,nb,j,nk,s,1)*jnred(nb,nk)*powneg1(-kb)
end
function sz_op(j::Real,s::Real,qns::Array{Int,2},p::Int)
   l = size(qns,1)
   out = spzeros(l,l)
   if s≠zero(s)&&p≠0
      for a ∈ 1:l, b ∈ a:l
         if abs(qns[a,1]-qns[b,1])≤1 && qns[a,2]==qns[b,2]
            @views out[b,a] = sqpart(j,s,0,qns[b,:],qns[a,:])
         end
      end
      dropzeros!(out)
      out .*= nred(s)*powneg1(s+j+1)
      out = Symmetric(out,:L)^p
   else
      out[diagind(out)] .+= 1.0
   end
   return out
end

function sq_op(j,s,q,qns)::SparseMatrixCSC{Float64, Int64}
   l = size(qns,1)
   out = spzeros(l,l)
   if s≠zero(s)#&&p≠0
      for a ∈ 1:l, b ∈ 1:l
         if abs(qns[a,1]-qns[b,1])≤1 && (q+qns[a,2]-qns[b,2])==0
            @views out[b,a] = sqpart(j,s,q,qns[b,:],qns[a,:])
         end
      end
      dropzeros!(out)
      out .*= nred(s)*powneg1(s+j+1+δ(1,q))*√2
      #the √2 is included to convert from spherical to cylinderical 
   else
      out[diagind(out)] .+= 1.0
   end
   return out
end

function sp_op(j::Real,s::Real,qns::Array{Int,2},p::Int
         )::SparseMatrixCSC{Float64, Int64} 
   if p≠0
      return sq_op(j,s,1,qns)^p
   else
      return spdiagm(ones(size(qns,1)))
   end
end
function sm_op(j::Real,s::Real,qns::Array{Int,2},p::Int
         )::SparseMatrixCSC{Float64, Int64} 
   return sq_op(j,s,-1,qns)^p
end
function spm_op(j::Real,s::Real,qns::Array{Int,2},p::Int
         )::SparseMatrixCSC{Float64, Int64} 
   return sp_op(j,s,qns,p) + sm_op(j,s,qns,p)
end

pa_op(ms::Array{Int},p::Int)::Diagonal{Float64, Vector{Float64}} = Diagonal(ms.^p)
function cos_op(ms::Array{Int},p::Int)::SparseMatrixCSC{Float64, Int64}
   if p==0
      out = I(size(ms,1))
   else
      out = fill(0.5, length(ms)-p)
      out = spdiagm(-p=>out, p=>out)
   end
   return out
end
function sin_op(ms::Array{Int},p::Int)::SparseMatrixCSC{Float64, Int64}
   #this is actually sin/2i as we are moving immediately multiplying it by
   #the i from Ny
   if p==0
      out = I(size(ms,1))
   else
      out = fill(0.25, length(ms)-p)
      out = spdiagm(-p=>out, p=>out)
   end
   return out
end

####   collected operators   ####
function rsr_op(j::Real,s::Real,qns::Array{Int,2},a::Int,b::Int,
         c::Int,d::Int,e::Int,f::Int,jp::Int)::SparseMatrixCSC{Float64, Int64}
   out = nnss_op(j,s,qns,a,b)*nz_op(qns,c)
   out = out*sz_op(j,s,qns,d)
   out = out*tplus!(np_op(qns,e)*sp_op(j,s,qns,f))
   out = out*iny_op(qns,jp)
   return dropzeros!(out)
end

tor_op(ms,g,h,j)::SparseMatrixCSC{Float64, Int64} = pa_op(ms,g)*
                                                      cos_op(ms,h)*sin_op(ms,j)

function tsr_op(prm::Float64,j::Real,s::Real,qns::Array{Int,2},ms::Array{Int},
                  a::Int,b::Int,c::Int,d::Int,e::Int,f::Int,g::Int,h::Int,
                  jp::Int)::SparseMatrixCSC{Float64, Int64}
   out = 0.25*prm*rsr_op(j,s,qns,a,b,c,d,e,f,jp)
   out = kron(tor_op(ms,g,h,jp),out)
   return dropzeros!(tplus!(out))
end
function tsr_op(prm::Float64,j::Real,s::Real,qns::Array{Int,2},
            ms::Array{Int},plist::Array{Int})::SparseMatrixCSC{Float64, Int64}
   #1/2 from a+a', 1/2 from np^0 + nm^0
   out = 0.25*prm*rsr_op(j,s,qns,plist[1],plist[2],plist[3],
                           plist[4],plist[5],plist[6],plist[9])
   out = kron(tor_op(ms,plist[7],plist[8],plist[9]),out)
   return dropzeros!(tplus!(out))
end

####   final construction and collection functions   ####

function hjbuild(sof,cdf::Array,cdo::Array,j,s,nf,tormat,ms)::SparseMatrixCSC{Float64, Int64}
   qns = qngen(j,s)
   #ms = msgen(nf,mc,σ)
   ℋ = hrot2(sof[1:4],qns) 
   if s==0.5
      ℋ .+= hsr(sof[5:9],j,s,qns)
   elseif s≥1
      ℋ .+= hsr(sof[5:9],j,s,qns)
      ℋ .+= hqu(sof[10:12],j,s,qns)
   end
   ℋ = kron(I(length(ms)), ℋ)
   ℋ .+= kron(tormat, I(size(qns,1)))
   if (s≥0.5)&&(nf≠0)
   ℋ .+= kron(pa_op(ms,1), sof[14]*nz_op(qns,1) + sof[15]*npm_op(qns,1) + 
               sof[17]*sz_op(j,s,qns,1) + sof[18]*spm_op(j,s,qns,1))
   elseif (s==zero(s))&&(nf≠0)
      ℋ += kron(pa_op(ms,1), sof[14]*nz_op(qns,1)+ sof[15]*npm_op(qns,1))
   else
   end
   for i in 1:length(cdf)
      ℋ .+= tsr_op(cdf[i],j,s,qns,ms,cdo[:,i] )
   end
   return dropzeros!(ℋ)
end

function tsrdiag(ctrl,sof,cdf,cdo,tormat,ms,nf,mcalc,j,s,σ,vtm)
   #println("ℋ time for J=$j")
   H = hjbuild(sof,cdf,cdo,j,s,nf,tormat,ms)
   if true ∈ isnan.(H)
      @warn "FUCK!!! j=$j, σ=$σ, NaN in H"
   end
   U = kron(ut(mcalc,σtype(nf,σ)),ur(j,s))
   H = (U*H*U)
   ### All lines commented with ### are for the Jacobi routine
   ###perm = kperm(j,s,mcalc)
   ###H = permute!(H,perm,perm)
   ###H, rvecs = jacobisparse(H, 3)#Int(j+s)+mcalc)
   ###rvecs = U*rvecs
   vals, vecs = eigen!(Symmetric(Matrix(H)))
   #@show vals[1:2*Int(2j+1)] ./csl
   ###perm = assignperm(vecs)
   if ctrl["assign"]=="RAM36"
      perm = ramassign(vecs,j,s,mcalc,σtype(nf,σ),vtm)
      vals = vals[perm]
      vecs = vecs[:,perm]
   elseif ctrl["assign"]=="expectk"
      vals, vecs = expectkassign!(vals,vecs,j,s,nf,mcalc,σ)      
   elseif ctrl["assign"]=="eeo"
      vals, vecs = eeoassign!(vals,vecs,j,s,nf,mcalc,σ)
   else
      vals, vecs = expectassign!(vals,vecs,j,s,nf,mcalc,σ)
   end
   ###vecs = rvecs*vecs 
   vecs = U*vecs
   return vals, vecs
end

function tsrcalc(ctrl,prm,stg,cdo,nf,vtm,mcalc,jlist,s,sd,σ)
   sof = prm[1:18]
   cdf = prmsetter(prm[19:end],stg)
   tormat, ms = htor2v2(sof[13:16],nf,mcalc,σ)
   msd = sd*Int(2*mcalc+(σtype(nf,σ)==2)+1)
   vtd = Int(vtm+1)
   σt = σtype(nf,σ)
   jmin = 0.5*iseven(sd)
   jmax = jlist[end]
   jfd = sd*Int(sum(2.0 .* collect(Float64,jmin:jmax) .+ 1.0))
   #msd = sd*mcd
   #mstrt, mstop = mslimit(nf,mcalc,σ)
   outvals = zeros(Float64,jfd*vtd)
   outpasz = zeros(Float64,jfd*vtd)
   outquns = zeros(Int,jfd*vtd,6)
   outvecs = zeros(Float64,Int((2*jmax+1)*msd),jfd*vtd)
   for j in jlist #thread removed for troubleshooting purposes
#   @threads for j in jlist
      jd = Int(2.0*j) + 1
      ###pull behavior should be move into TSRDIAG moving forward
      ###pull = indpuller(vtm,mcalc,σt,Int(jd*sd))
      sind, find = jvdest(j,s,vtm) 
      tvals, tvecs = tsrdiag(ctrl,sof,cdf,cdo,tormat,ms,nf,mcalc,j,s,σ,vtm)
      outvals[sind:find] = tvals
      outquns[sind:find,:] = qnlabv(j,s,nf,vtm,σ)
      outvecs[1:jd*msd,sind:find] = tvecs###[:,pull]
   end
   return outvals, outvecs, outquns
end

function tsrcalc2(prm,stg,cdo,nf,ctrl,jlist)
   s = ctrl["S"]
   mcalc = ctrl["mcalc"]
   vtm = ctrl["vtmax"]
   sof = prm[1:18]
   cdf = prmsetter(prm[19:end],stg)
   vtd = Int(vtm+1)
   jmin = 0.5*iseven((2*s+1))
   jmax = jlist[end,1]
   jfd = Int((2s+1)*sum(2. .*collect(Float64,jmin:jmax).+1.))
   σcnt = σcount(nf)
   fvls = zeros(Float64,jfd*vtd,σcnt)
   fqns = zeros(Int,jfd*vtd,6,σcnt)
   fvcs = zeros(Float64,Int((2*s+1)*(2*jmax+2)*(2*mcalc+(iseven(nf))+1)),jfd*vtd,σcnt)
   @time for sc in 1:σcnt
      σ = sc - 1
      tormat, ms = htor2v2(sof[13:16],nf,mcalc,σ)
      msd = Int((2*mcalc+(σtype(nf,σ)==2)+1)*(2s+1))
      jmsd = Int(msd*(2*jmax+1))
      jsvd = Int(jfd*vtd)
      jsublist = jlist[isequal.(jlist[:,2],σ), 1] .* 0.5
      @threads for j in jsublist
         jd = Int(2.0*j) + 1
         sind, find = jvdest(j,s,vtm) 
         fvls[sind:find,sc], fvcs[1:jd*msd,sind:find,sc] = tsrdiag(ctrl,
            sof,cdf,cdo,tormat,ms,nf,mcalc,j,s,σ,vtm)
         fqns[sind:find,:,sc] = qnlabv(j,s,nf,vtm,σ)
      end
   end
   return fvls, fvcs, fqns
end
