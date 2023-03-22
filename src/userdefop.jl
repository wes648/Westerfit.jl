
#=
include("@__DIR__/WIGXJPF.jl")
include("@__DIR__/jacobi.jl")
include("@__DIR__/filehandling.jl")
using LinearAlgebra, SparseArrays, Base.Threads, .WIGXJPF, LinearAlgebra.BLAS
include("@__DIR__./userdefop.jl")
=#

indmax(vtm,σt)::Int = ceil(vtm/2)*cospi(vtm) + (vtm≤zero(vtm))*δi(σt,2)
function vtlist(vtm,σt) 
   if vtm == zero(vtm)
      return [indmax(vtm,σt)]
   else
      return sort([indmax(vtm,σt); indmax(vtm-1,σt)])
   end
end
function vtcoll(vtm,σt) 
   if vtm == zero(vtm)
      return [indmax(vtm,σt)]
   else
      a = indmax(vtm,σt)
      b = indmax(vtm-1,σt)
      ref = sort([a; b])
      return collect(ref[1]:ref[2])
   end
end
function indpuller(vtm,mc,σt,jsd)
   a = indmax(vtm,σt)
   b = indmax(vtm-1,σt)
   ref = sort([a; b])
   if σt != 2
      c = jsd*(mc + ref[1]) + 1
      d = jsd*(mc + ref[2] +1)
   else
      c = jsd*(mc + ref[1] - (vtm>zero(vtm))) + 1
      d = jsd*(mc + ref[2] + 1 - (vtm>zero(vtm)))
   end
   return collect(c:d)
end
function jvdest(j,s,vtm)
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return snd,fnd
end
function jvlinds(j,s,vtm)
   snd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return collect(snd:fnd)
end
function jinds(j,s,m,σt)
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (2*m+1+δ(σt,2))*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (2*m+1+δ(σt,2))*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return snd,fnd
end
function jlinds(j,s,m,σt)
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (2*m+1+δ(σt,2))*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (2*m+1+δ(σt,2))*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return collect(snd:fnd)
end

function Δlist(J,S)
   max = Int(J+S)
   min = Int(abs(J-S))
   return collect(min:max)
end

function kperm(n::Int)::Array{Int}
   sortperm(Int.(cospi.(collect(-n:n).+isodd(n))) .* collect(-n:n))
end
function kperm(j,s)::Array{Int}
   perm = zeros(Int,Int((2*j+1)*(2*s+1)))
   shift = 0
   for n in Δlist(j,s)
      nd = 2*n+1
      perm[(1+shift):(nd+shift)] = kperm(n) .+ shift
      shift += nd
   end
   return perm
end
function kperm(j,s,shift::Int,jsd::Int,Δl::Array,perm)::Array{Int}
   #perm = zeros(Int,jsd)
   for n in Δl
      nd = 2*n+1
      perm[(1+shift):(nd+shift)] = kperm(n) .+ shift
      shift += nd
   end
   return perm#, shift
end
function kperm(j,s,m)
   jsd = Int((2*j+1)*(2*s+1))
   Δlst = Δlist(j,s)
   shift = 0
   perm = zeros(Int,jsd*(2*m+1))
   for i in 0:(2*m)
      perm = kperm(j,s,shift,jsd,Δlst,perm)
      shift += jsd
   end
   return perm
end


function qngen(n,nf,m,σ)
   nd = Int(2*n+1)
   md = Int(2*m+(σtype(nf,σ)==2)+1)
   narray = fill(n,nd*md)
   karray = kron(ones(Int,md),collect(Int,-n:n))
   kcarray = k2kc.(narray,karray)
   marray = kron(msbuilder(nf,m,σ),ones(Int,nd))
   σarray = fill(σ,nd*md)
   out = hcat(narray,karray,kcarray,marray,σarray)
end
function qngen(j,s,nf,m,σ)
   nlist = Δlist(j,s)
   jsd = Int((2*j+1)*(2*s+1))
   md = Int(2*m+(σtype(nf,σ)==2)+1)
   out = zeros(Int,0,3)
   for n in nlist
      nd = Int(2*n+1)
      part = zeros(Int,nd,3)
      part[:,1] = fill(n,nd)
      part[:,2] = collect(Int,-n:n)
      part[:,3] = k2kc.(part[:,1],part[:,2])
      out = vcat(out,part)
   end
   out = kron(ones(Int,md),out)
   marray = kron(msbuilder(nf,m,σ),ones(Int,jsd))
   out = hcat(fill(Int(2*j),size(out,1)),out,marray,fill(σ,jsd*md))
   return out
end
function qngenv(j,s,nf,vtm,σ)
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
   out = kron(ones(Int,vd),out)
   vtrray = kron(nf .* vtcoll(vtm,σt) .+ σ,ones(Int,jsd))
   out = hcat(fill(Int(2*j),size(out,1)),out,vtrray,fill(σ,jsd*vd))
   return out
end
function k2kc(n,k)
   ka = abs(k)
   if k < 0
      kc = n - ka + 1
   else
      kc = n - ka
   end
   return kc
end

eh(x::Int)::Float64 = √(x*(x+1))
eh(x::Float64)::Float64 = √(x*(x+1.0))
assignperm(vec) = sortperm([iamax(vec[:,i]) for i in 1:size(vec,2)])

function fh(x::Int,y::Int)::Float64
   out = □rt((x-y)*(x+y+1))
   return out
end
function fh(x::Float64,y::Float64)::Float64
   out = □rt((x-y)*(x+y+1))
   return out
end
function gh(x,y)::Float64
   out = □rt((x-y)*(x-y-1.0))
   return out
end
function □rt(x)::Float64
   if x ≤ zero(x)
      return 0.0
   else
      return √x
   end
end
function θ(j,n,s)
   if s==zero(s)
      out = 0.0
   elseif s==0.5
      out = (n-j)/(j+0.5)
   else
      out = n*(n+1.0) + s*(s+1.0) - j*(j+1.0)
      if out != zero(out)
         out = out/(2.0*n*(n+1.0))
      end
   end
   return out
end
function ϕ(j,n,s)
   if s==zero(s)
      out = 0.0
   elseif s==0.5
      out = -1.0/(j+0.5)
   else
      out = (n-j+s)*(n+j+s+1)*(s+j-n+1)*(n+j-s)
      out *= 1.0/((2.0*n-1.0)*(2.0*n+1.0))
      out = -sqrt(out)/n
   end
   return out
end

δ(x::Int,y::Int)::Float64 = x==y
δ(x::Float64,y::Float64)::Float64 = x==y
δ(x::Float64,y::Int)::Float64 = x==convert(Float64,y)
δ(x::Int,y::Float64)::Float64 = convert(Float64,x)==y
δi(x::Int,y::Int)::Int = x==y

T(l::Int,q::Int)::Int = l*(l+1) + q + 1
Tμ(q::Int)::Int = q + 2
Tq(q::Int)::Int = abs(q) + 1 #quadrupole variant (only has 2nd rank components)
Tsr(l::Int,q::Int)::Int = δi(l,2) + abs(q) + 1 #Cs sr version, no 1st rk, & symm

function kgen(n::Int)::Array{Int,2}
   return kron(collect(-n:n),ones(Int,1,2*n+1))
end
function kgen2(n::Int)::LowerTriangular{Int, Matrix{Int}}
   ks = collect(-n:n)
   out = LowerTriangular(zeros(Int,length(ks),length(ks)))
   @simd for i in 1:length(ks)
      @inbounds out[i:end,i] = ks[i:end]
   end
   return out
end

function kgeni(n::Int,lb::Int)::Array{Int,2}
   return kron(collect(-n:n),ones(Int,1,lb))
end
function kgeni(j::Float64,s::Float64,lb::Int)::Array{Int,2}
   ns = Δlist(j,s)
   out = zeros(Int,Int((2.0*s+1.0)*(2.0*j+1.0)))
   si = 1
   for n in ns 
      fi = si + 2*n
      out[si:fi] = collect(-n:n)
      si = fi + 1
   end
   return kron(out,ones(Int,1,lb))
end
function kgen(j::Float64,s::Float64)::Array{Int,2}
   ns = Δlist(j,s)
   out = zeros(Int,Int((2.0*s+1.0)*(2.0*j+1.0)))
   si = 1
   for n in ns 
      fi = si + 2*n
      out[si:fi] = collect(-n:n)
      si = fi + 1
   end
   lb = convert(Int,(2.0*j+1.0)*(2.0*s+1.0))
   return kron(out,ones(Int,1,lb))
end
function ngen(n::Int)::Array{Int,2}
   return kron(fill(n,2*n+1),ones(Int,1,2*n+1))
end
function ngeni(j::Float64,s::Float64,lb::Int)::Array{Int,2}
   ns = Δlist(j,s)
   out = zeros(Int,Int((2.0*s+1.0)*(2.0*j+1.0)))
   si = 1
   for n in ns 
      fi = si + 2*n
      out[si:fi] = fill(n,2*n+1)
      si = fi + 1
   end
   return kron(out,ones(Int,1,lb))
end
function ngen(j::Float64,s::Float64)::Array{Int,2}
   ns = Δlist(j,s)
   lb = convert(Int, (2.0*s+1.0)*(2.0*j+1.0))
   out = zeros(Int,lb)
   si = 1
   for n in ns 
      fi = si + 2*n
      out[si:fi] = fill(n,2*n+1)
      si = fi + 1
   end
   return kron(out,ones(Int,1,lb))
end

σcount(nfold::Int)::Int = floor(Int,nfold/2)+1
function σtype(nfold,σ)
   if σ==zero(σ) # A state
      return 0
   elseif (iseven(nfold))&&(σ==(σcount(nfold)-1)) # B state
      return 2
   else # E state
      return 1
   end
end
function msbuilder(T::Type,nfold::Number,mcalc::Number,σ::Number)
   σt = σtype(nfold,σ)
   lim = mcalc*nfold
   if σt==0
      marray = collect(T,-lim:nfold:lim)
   elseif σt==2
      lim += σ
      marray = collect(T,-lim:nfold:lim)
   else
      marray = collect(T,(-lim+σ):nfold:(lim+σ))
   end
   return marray
end
msbuilder(nfold::Int,mcalc::Int,σ::Int)::Array{Int} = msbuilder(Int,nfold,mcalc,σ)
mgen(nf::Int,mc::Int,σ::Int)::Array{Int,2} = kron(msbuilder(nf,mc,σ), ones(Int,1,2*mc+1))

function mslimit(nfold,mcalc,σ)::Tuple{Int, Int}
   σt = σtype(nfold,σ)
   lim = mcalc*nfold
   if σt==0
      return -lim, lim
   elseif σt==2
      lim += σ
      return -lim, lim
   else
      return (-lim+σ), (lim+σ)
   end
end

nred(n::Int)::Float64 = √(n*(n+1)*(2*n+1))
nred(n::Float64)::Float64 = √(n*(n+1.0)*(2.0*n+1.0))
nnred(n::Int)::Float64 = n*(n+1)*(2*n+1)
nnred(n::Float64)::Float64 = n*(n+1.0)*(2.0*n+1.0)

function hrred(n::Int,l::Int)::Float64
   return wig6j(n,n,1,l,1,n)*nnred(n)*√(2*l+1)*(-1)^l
end
function hrelem(pr::Float64,l::Int,q::Int,nb,kb,nk,kk)::Array{Float64,2}
   @. return pr*wig3j(nb,l,nk,-kb,q,kk)*hrred(nk,l)*δ(nb,nk)*(-1)^(nk-kb)
end
function hrdiag(bn,bk,nb,kb,nk,kk)::Array{Float64}
   @. return (bn*nk*(nk+1.0) + bk*kk^2)*δ(nb,nk)*δ(kb,kk)
end
function hrfu1(dab,nb,kb,nk,kk)::Array{Float64}
   @. return dab*(kk+0.5)*fh(nk,kk)*δ(kb,kk+1)*δ(nb,nk)
end
function hrfl1(dab,nb,kb,nk,kk)::Array{Float64}
   @. return dab*(kk-0.5)*fh(nk,kk-1)*δ(kb,kk-1)*δ(nb,nk)
end
function hrfu2(bpm,nb,kb,nk,kk)::Array{Float64}
   @. return bpm*fh(nk,kk)*fh(nk,kk+1)*δ(kb,kk+2)*δ(nb,nk)
end
function hrfl2(bpm,nb,kb,nk,kk)::Array{Float64}
   @. return bpm*fh(nk,kk-1)*fh(nk,kk-2)*δ(kb,kk-2)*δ(nb,nk)
end

function hrdiag(bn::Float64,bk::Float64,nk,kk)::Array{Float64}
   @. return bn*nk*(nk+1.0) + bk*kk^2
end
function hrfu1(dab::Float64,nk,kk)::Array{Float64}
   @. return dab*(kk+0.5)*fh(nk,kk)
end
function hrfu2(bpm::Float64,nk,kk)::Array{Float64}
   @. return bpm*fh(nk,kk)*fh(nk,kk+1)
end

function hrot2s(pr,nb,kb,nk,kk)
   out = spzeros(size(kb))
   lst = diagind(out)
   out[lst] .= hrdiag(pr[2],pr[1],nb[lst],kb[lst],nk[lst],kk[lst])
   lst = diagind(out,1)
   out[lst] .= hrfu1(pr[4],nb[lst],kb[lst],nk[lst],kk[lst])
   lst = diagind(out,-1)
   out[lst] .= hrfl1(pr[4],nb[lst],kb[lst],nk[lst],kk[lst])
   lst = diagind(out,2)
   out[lst] .= hrfu2(pr[3],nb[lst],kb[lst],nk[lst],kk[lst])
   lst = diagind(out,-2)
   out[lst] .= hrfl2(pr[3],nb[lst],kb[lst],nk[lst],kk[lst])
   return out
end
function hrot2v2(pr,nb,kb,nk,kk)::SparseMatrixCSC{Float64, Int64}#fastest
   out = spzeros(Float64, size(kb))
   lst = diagind(out,1)
   out[lst] .= hrfu1(pr[4],nk[lst],kk[lst])
   lst = diagind(out,2)
   out[lst] .= hrfu2(pr[3],nk[lst],kk[lst])
   out .+= transpose(out)
   lst = diagind(out)
   out[lst] .= hrdiag(pr[2],pr[1],nk[lst],kk[lst])
   return out
end
function hrtest(n)::SparseMatrixCSC{Float64, Int64}
   pr = [1.75; 1.25; 0.25; 0.02]
   nk = ngen(n)
   kk = kgen(n)
   nb = Matrix(transpose(nk))
   kb = Matrix(transpose(kk))
   return hrot2v2(pr,nb,kb,nk,kk)
end

function hrlpart(out,pr,l::Int,nb,kb,nk,kk)::Array{Float64,2}
   @simd for q in -l:l
      out .+= hrelem(pr[T(l,q)],l,q,nb,kb,nk,kk)
   end
   return out
end
function hrot(pr,n)
   nks = ngen(n)
   nbs = Matrix(transpose(nks))
   kks = kgen(n)
   kbs = transpose(kks)
   out = zeros(size(kks))
   @simd for l in 0:2:2
      out .= hrlpart(out,pr,l,nbs,kbs,nks,kks)
   end
   return out
end
function hrot3(pr,nb,kb,nk,kk)
   out = zeros(size(kk))
   @simd for l in 0:2:2
      out .= hrlpart(out,pr,l,nb,kb,nk,kk)
   end
   return out
end
function hrot(pr,j,s)
   lb = convert(Int,(2.0*s+1.0)*(2.0*j+1.0))
   nks = ngen(j,s,lb)
   nbs = Matrix(transpose(nks))
   kks = kgen(j,s,lb)
   kbs = transpose(kks)
   out = zeros(size(kks))
   @simd for l in 0:2:2
      out .= hrlpart(out,pr,l,nbs,kbs,nks,kks)
   end
   return out
end

function nsred(l::Int,nb,nk)
   @. return 0.5*(nred(nk)*wig6j(1,1,l,nb,nk,nk) + nred(nb)*wig6j(1,1,l,nk,nb,nb))
end
function jsred(j,s,nb,nk)
   @. return wig6j(nk,s,j,s,nb,1)*jnred(nb,nk)*nnred(s)
end
function srelem(pr::Float64,l::Int,q::Int,j,s,nb,kb,nk,kk)#::Array{Float64,2}
   @. return pr*wig3j(nb,l,nk,-kb,q,kk)*√(2.0*l+1.0)*
             nsred(l,nb,nk)*jsred(j,s,nb,nk)*(-1.0)^(j+s-kb+δ(-1,q))
end
function srlpart(out,pr,l::Int,j,s,nb,kb,nk,kk)#::Array{Float64,2}
   @simd for q in -l:l
      out += srelem(pr[Tsr(l,q)],l,q,j,s,nb,kb,nk,kk)
   end
   return out
end
function hsr(pr,j,s)
   lb = convert(Int,(2.0*s+1.0)*(2.0*j+1.0))
   nks = ngen(j,s)
   nbs = Matrix(transpose(nks))
   kks = kgen(j,s)
   kbs = transpose(kks)
   out = zeros(size(kks))
   @simd for l in 0:2:2
      out .= srlpart(out,pr,l,j,s,nbs,kbs,nks,kks)
   end
   return out
end
function hrsr(rpr,spr,j,s)
   lb = convert(Int,(2.0*s+1.0)*(2.0*j+1.0))
   nks = ngen(j,s,lb)
   nbs = Matrix(transpose(nks))
   kks = kgen(j,s,lb)
   kbs = transpose(kks)
   out = zeros(size(kks))
   @simd for l in 0:2:2
      out .= hrlpart(out,rpr,l,nbs,kbs,nks,kks)
      out .= srlpart(out,spr,l,j,s,nbs,kbs,nks,kks)
   end
   return out
end
function hrsr(rpr,spr,qpr,j,s,nb,kb,nk,kk)
   out = hrot2v2(rpr,nb,kb,nk,kk)
   if s == zero(s)
   elseif s < one(s)
      @simd for l in 0:2:2
         out .= srlpart(out,spr,l,j,s,nb,kb,nk,kk)
      end
   else
      @simd for l in 0:2:2
         out .= srlpart(out,spr,l,j,s,nb,kb,nk,kk)
         out .= qulpart(out,qpr,j,s,nb,kb,nk,kk)
      end
   end
   return out
end
function hsrq(out,spr,qpr,j,s,nb,kb,nk,kk)
   @simd for l in 0:2:2
      out .= srlpart(out,spr,l,j,s,nb,kb,nk,kk)
      out .= qulpart(out,qpr,j,s,nb,kb,nk,kk)
   end
   return out
end

function wigdiv(x::Array,s::Number)::Array
   if s<one(s)
      return zero(x)
   else
      return x ./ wig3j(s,2,s,-s,0,s)
   end
end
function qured(j,s,nb,nk)
   @. return 0.25*jnred(nb,nk)*wig6j(j,s,nb,2,nk,s)
end
function quelem(pr,q,j,s,nb,kb,nk,kk)#::Array{Float64,2}
   @. return pr*qured(j,s,nb,nk)*
             wig3j(nb,2,nk,-kb,q,kk)*(-1.0)^(nk+nb-kk+s+j+1+δ(-1,q))
end
function qulpart(out,pr,j,s,nb,kb,nk,kk)#::Array{Float64,2}
   @simd for q in -2:2
      out += quelem(pr[Tq(q)],q,j,s,nb,kb,nk,kk)
   end
   out = wigdiv(out,s)
   return out
end

function jnred(j::Float64,n::Int)::Float64
   return √((2.0*j+1.0)*(2*n+1))
end
function jnred(j::Int,n::Int)::Float64
   return √((2*j+1)*(2*n+1))
end

function cart2sphr(inp::Array{Float64,2})::Array{Float64,1}
   out = zeros(9)
   out[1] = -sum(diag(inp))/√3
   out[2] = 0.5*(inp[1,2] - inp[2,1])
   out[3] = 0.0
   out[4] = 0.5*(inp[1,2] - inp[2,1])
   out[5] = 0.5*(inp[2,2] - inp[3,3])
   out[6] = 0.5*(inp[1,2] + inp[2,1])
   out[7] = (3.0*inp[1,1] - sum(diag(inp)))/√6
   out[8] = -0.5*(inp[1,2] + inp[2,1])
   out[9] = 0.5*(inp[2,2] - inp[3,3])
   return out
end
function cart2sphr(inp::Array{Float64,1})::Array{Float64,1}
   out = zeros(3)
   out[1] = (inp[3]+inp[1])/√2
   out[2] = inp[2]
   out[3] = -(inp[3]-inp[1])/√2
   return out
end

#function nzop(p::Int,nb::Array{Int,2},kb::Array{Int,2},
#              nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
#   return @. δ(nb,nk)*δ(kb,kk)*kk^p
#end
function nzop(p::Int,nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::Diagonal{Float64, Vector{Float64}}
   return Diagonal(kk) ^ p
end
function ntop(p::Int,nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::Diagonal{Float64, Vector{Float64}}
   return eh.(Diagonal(nk)) ^ p
end
function nyel(nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}#::Array{Float64,2}
   return @. δ(nb,nk)*(δ(kb-1,kk)*fh(nk,kk) - δ(kb+1,kk)*fh(nk,kk-1))
end
function nyop(p::Int,nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}#::Array{Float64,2}
   return nyel(nb,kb,nk,kk)^p
end
function nmel(nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}
   return @. δ(nb,nk)*δ(kb,kk+1)*fh(nk,kk)
end
function npel(nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}
   return @. δ(nb,nk)*δ(kb,kk-1)*fh(nk,kk-1)
end
function npmp(p::Int,nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}
   out = npel(nb,kb,nk,kk)^p
   return out + transpose(out)
   #return npel(nb,kb,nk,kk)^p + nmel(nb,kb,nk,kk)^p
end
function rotop(pr::Float64,t::Int,z::Int,x::Int,
               nb::Array{Int,2},kb::Array{Int,2},
               nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   tz = ntop(t,nb,kb,nk,kk).*nzop(z,nb,kb,nk,kk) #./ 2.0
   pm = npmp(x,nb,kb,nk,kk)
   return (pm*tz + tz*pm) .* (pr/4.0)
end
function nsop(p::Int,j::Float64,s::Float64,nb::Array{Int,2},kb::Array{Int,2},
                     nk::Array{Int,2},kk::Array{Int,2})::Diagonal{Float64, Vector{Float64}}
   return (0.5 *((eh(j)-eh(s))*I(size(nk,1)) - eh.(Diagonal(nk)))) ^ p
end
function nsop_old(p::Int,j::Float64,s::Float64,nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   return @. (0.5*δ(nb,nk)*δ(kb,kk)*(eh(j) - eh(nk) - eh(s)))^p
end
function szen(j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})
   return @. δ(nb,nk)*δ(kb,kk)*kk*θ(j,nk,s)
end
function szem(j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})
   return @. δ(nb+1,nk)*δ(kb,kk)*□rt(nk^2 - kk^2)*ϕ(j,nk,s)*0.5
end
function szep(j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})
   return @. δ(nb-1,nk)*δ(kb,kk)*□rt(nb^2 - kk^2)*ϕ(j,nb,s)*0.5
end
function szop(p::Int,j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})#::Array{Float64,2}
   if s==zero(s)
      return eye(size(nb,1))
   else
      out = szep(j,s,nb,kb,nk,kk)
      out += transpose(out)
      out += szen(j,s,nb,kb,nk,kk)
   return out
   end
end

#These functions, those more beautiful, don't work
function szelem(j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   return @. jnred(nb,nk)*wig6j(s,nb,j,nk,s,1)*wig3j(nb,1,nk,-kb,0,kk)*nred(s)*(-1)^(s+j-kb+1)
end
function szop2(p::Int,j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   return szelem(j,s,nb,kb,nk,kk)^p
end
function szelem3(j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}
   return @. jnred(nb,nk)*wig6j(s,nb,j,nk,s,1)*wig3j(nb,1,nk,-kb,0,kk)*nred(s)*(-1)^(s+j-kb+1)
end
function szop3(p::Int,j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}
   return szelem3(j,s,nb,kb,nk,kk)^p
end


#function rsrop(pr::Float64,a::Int,b::Int,c::Int,d::Int,e::Int,
#               j,s,nb::Array{Int,2},kb::Array{Int,2},
#               nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
#   tz = ntop(a,nb,kb,nk,kk).*nzop(b,nb,kb,nk,kk).*nsop(d,j,s,nb,kb,nk,kk) #these all commute!
#   pm = npmp(c,nb,kb,nk,kk)
#   sz = szop(e,j,s,nb,kb,nk,kk)
#   return (tz*sz*pm + pm*sz*tz) .* (pr/4.0)
#end
function rsrop(a::Int,b::Int,c::Int,d::Int,e::Int,h::Int,
               j,s,nb::Array{Int,2},kb::Array{Int,2},
               nk::Array{Int,2},kk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}
   #::SparseMatrixCSC{Float64, Int64}#::Array{Float64,2}
   op = ntop(a,nb,kb,nk,kk)*nzop(b,nb,kb,nk,kk)*nsop(d,j,s,nb,kb,nk,kk) 
   #the above all commute!
   op *= sparse(npmp(c,nb,kb,nk,kk))
   op *= sparse(nyop(1-δi(0,h),nb,kb,nk,kk))
   op *= sparse(szop(e,j,s,nb,kb,nk,kk))
   return 0.25 * (op + transpose(op))
end

function paop(p::Int,mb::Array{Int,2},
             mk::Array{Int,2})::Diagonal{Float64, Vector{Float64}}
   return Diagonal(mk) ^ p
end
function cosp(p::Int,mb::Array{Int,2},mk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}
   return @. (δ(mb+p,mk)+δ(mb-p,mk))/2.0
end
function sinp(p::Int,mb::Array{Int,2},mk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}
   return @. (δ(mb+p,mk) + δ(mb-p,mk)) / 2.0
end
function torop(pr::Float64,p::Int,c::Int,s::Int,
               mb::Array{Int,2},mk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}
   op = paop(p,mb,mk)
   op *= sparse(cosp(c,mb,mk)*sinp(s,mb,mk))
   op .+= transpose(op)
   #out = symm('L','U',pa,sc) + symm('R','U',pa,sc)
   #out = (pa*sc + sc*pa)
   return (pr*0.5) * op
end
function torop(pr::Tuple,p::Tuple,c::Tuple,s::Tuple,
               mb::Array{Int,2},mk::Array{Int,2})#::Array{Float64,2}
   out = torop(pr[1],p[1],c[1],s[1],mb,mk)
   for i in 2:length(pr)
      out .+= torop(pr[1]*pr[i],p[i],c[i],s[i],mb,mk)
   end
   return out
end
function torop(pr::Float64,p::Int,c::Int,
               mb::Array{Int,2},mk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}
   op = paop(p,mb,mk)*sparse(cosp(c,mb,mk))
   return (op + transpose(op)) .* (pr/2.0)
end

function tsrop(pr::Float64,a::Int,b::Int,c::Int,d::Int,e::Int,f::Int,g::Int,h::Int,
               j::Float64,s::Float64,nb::Array{Int,2},kb::Array{Int,2},
               mb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2},
               mk::Array{Int,2})::SparseMatrixCSC{Float64, Int64}#::Array{Float64,2}
   return dropzeros!(kron(torop(pr,f,g,h,mb,mk),rsrop(a,b,c,d,e,h,j,s,nb,kb,nk,kk)))
end
function tsrop(pr::Tuple,a::Tuple,b::Tuple,c::Tuple,d::Tuple,e::Tuple,f::Tuple,
               g::Tuple,h::Tuple,
               j::Float64,s::Float64,nb::Array{Int,2},kb::Array{Int,2},
               mb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2},
               mk::Array{Int,2})#::Array{Float64,2}
   out = tsrop(pr[1],a[1],b[1],c[1],d[1],e[1],f[1],g[1],h[1],j,s,nb,kb,mb,nk,kk,mk)
   for i in 2:length(pr)
      out += tsrop(pr[1]*pr[i],a[i],b[i],c[i],d[i],e[i],f[i],g[i],h[i],
                   j,s,nb,kb,mb,nk,kk,mk)
   end
   return dropzeros!(out)
end
function tsrop(pr::Float64,op::Array,j::Float64,s::Float64,
               nb::Array{Int,2},kb::Array{Int,2},
               mb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2},
               mk::Array{Int,2})#::Array{Float64,2}
   return tsrop(pr,op[1],op[2],op[3],op[4],op[5],op[6],op[7],op[8],
                j,s,nb,kb,mb,nk,kk,mk)
end

function hbuild(rot,spi,qua,tor,cdf,cdo,j,s,mb,mk)
   nk = ngen(j,s)
   kk = kgen(j,s)
   nb = Matrix(transpose(nk))
   kb = Matrix(transpose(kk))
   hout = kron(eye(size(mk,1)), hrsr(rot,spi,qua,j,s,nb,kb,nk,kk))
   htor = torop(tor[1],2,0,0,mb,mk)
   htor .+= torop(tor[3],0,1,0,mb,mk).+ torop(tor[3],0,0,0,mb,mk)
   hout .+= kron(htor,eye(size(nk,1))) .+ 
                 tsrop(tor[2],0,1,0,0,0,1,0,0,j,s,nb,kb,mb,nk,kk,mk)
   hout .+= tsrop(tor[4],0,0,0,0,1,1,0,0,j,s,nf,nb,kb,mb,nk,kk,mk)
   @simd for i in 1:length(cdf)
      hout .+= tsrop(cdf[i],cdo[:,i],j,s,nb,kb,mb,nk,kk,mk)
   end
   return hout
end

function htor(sof,cdf::Array,cdo::Array,nf,mcalc,σ)
   mk = mgen(nf,mcalc,σ)
   mb = Matrix(transpose(mk))
   out = torop(sof[12],2,0,mb,mk)
   out += sof[14]*eye(size(out,1)) - torop(sof[14],0,1,mb,mk)
   @simd for i in 1:length(cdf)
      out += torop(cdf[i],cdo[6,i],cdo[7,i],cdo[8,i],mb,mk)
   end
   return out, mk, mb
end
function htor(sof,cdf::Nothing,cdo::Nothing,nf,mcalc,σ)
   mk = mgen(nf,mcalc,σ)
   mb = Matrix(transpose(mk))
   out = torop(sof[12],2,0,mb,mk)
   out += sof[14]*eye(size(out,1)) - torop(sof[14],0,1,mb,mk)
   return out, mk, mb
end
function htor(sof,nf,mcalc,σ) #this is the current one being called by the code
   mk = mgen(nf,mcalc,σ)
   mb = Matrix(transpose(mk))
   out = torop(sof[12],2,0,mb,mk)
   out += sof[14]*eye(size(out,1)) - torop(sof[14],0,nf,mb,mk)
   return out, mk, mb
end
function htorq(sof,nf,mb,mk)
   out = torop(sof[1],2,0,mb,mk)
   out += sof[3]*eye(size(out,1)) - torop(sof[3],0,nf,mb,mk)
   return out
end
function hjbuild(sof,cdf::Array,cdo::Array,tormat,j,s,mb,mk)
   nk = ngen(j,s)
   kk = kgen(j,s)
   nb = Matrix(transpose(nk))
   kb = Matrix(transpose(kk))
   #scale up tormat & add -2ρF
   hout = kron(tormat,eye(size(nk,1))) 
   hout -= kron(sof[13] * Diagonal(mk), Diagonal(kk))
   #tsrop(sof[13],0,1,0,0,0,1,0,0,j,s,nb,kb,mb,nk,kk,mk)
   #add 2nd order ro, spi, qua
   hout += kron(eye(size(mk,1)), hrsr(sof[1:4],sof[5:8],sof[9:11],j,s,nb,kb,nk,kk))
   #println("hout type = $(typeof(hout))")
   #if s != zero(s)#add η
   #hout += tsrop(sof[15],0,0,0,0,1,1,0,0,j,s,nb,kb,mb,nk,kk,mk)
   #end
   @simd for i in 1:length(cdf)
      hout += tsrop(cdf[i],cdo[:,i],j,s,nb,kb,mb,nk,kk,mk)
   end
   return hout
end
#=function hjbuild(sof,cdf::Nothing,cdo::Nothing,tormat,j,s,mb,mk)
   nk = ngen(j,s)
   kk = kgen(j,s)
   nb = Matrix(transpose(nk))
   kb = Matrix(transpose(kk))
   #scale up tormat & add -2ρF
   hout = kron(tormat,eye(size(nk,1))) -
               tsrop(sof[13],0,1,0,0,0,1,0,0,j,s,nb,kb,mb,nk,kk,mk)
   #add 2nd order ro, spi, qua
   hout += kron(eye(size(mk,1)), hrsr(sof[1:4],sof[5:8],sof[9:11],j,s,nb,kb,nk,kk))
   #add η
   hout += tsrop(sof[15],0,0,0,0,1,1,0,0,j,s,nb,kb,mb,nk,kk,mk)
   return hout
end=#

function tsrdiag(sof,cdf,cdo,tormat,nf,mcalc,mb,mk,j,s,σ,σt)
   #fuse sof & cdf in tsdriag call 
   H = hjbuild(sof,cdf,cdo,tormat,j,s,mb,mk)
   #println("H type = $(typeof(H))")
   if σtype(nf,σ) != 1
      U = ur(j,s,mcalc,σt)*ut(mcalc,σt,j,s)
   else
      U = ur(j,s,mcalc,σt)
   end
   H = (U*H*U)
   perm = kperm(j,s,mcalc)
   H = permute!(H,perm,perm)
   #H, rvecs = SparseSweep(H)
   H, rvecs = limsparsweep(H,3)
   rvecs = U*rvecs
   vals, vecs = LAPACK.syev!('V', 'U', Matrix(H))
   perm = assignperm(vecs)
   vals = vals[perm]
   vecs = vecs[perm,:]
   vecs = rvecs*vecs 
   return vals, vecs
end

function tsrcalc(prm,stg,cdo,nf,vtm,mcalc,jlist,s,sd,σ)
   sof = prm[1:15]
   cdf = prmsetter(prm[16:end],stg)
   tormat, mk, mb = htor(sof,nf,mcalc,σ)
   #sd = Int(2.0*s+1.0)
   mcd = Int(2*mcalc+(σtype(nf,σ)==2)+1)
   vtd = Int(vtm+1)
   σt = σtype(nf,σ)
   jmin = 0.5*iseven(sd)
   jmax = jlist[end]
   jfd = sd*Int(sum(2.0 .* collect(Float64,jmin:jmax) .+ 1.0))
   msd = sd*mcd
   mstrt, mstop = mslimit(nf,mcalc,σ)
   outvals = zeros(Float64,jfd*vtd)
   outquns = zeros(Int,jfd*vtd,6)
   outvecs = zeros(Float64,Int(sd*(2*jmax+1)*mcd),jfd*vtd)
#   for j in jlist #thread removed for troubleshooting purposes
   @threads for j in jlist
      jd = Int(2.0*j) + 1
      pull = indpuller(vtm,mcalc,σt,Int(jd*sd))
      sind, find = jvdest(j,s,vtm) 
      tvals, tvecs = tsrdiag(sof,cdf,cdo,tormat,nf,mcalc,mb,mk,j,s,σ,σt)
      outvals[sind:find] = tvals[pull]
      outquns[sind:find,:] = qngenv(j,s,nf,vtm,σ)
      outvecs[1:jd*msd,sind:find] = tvecs[:,pull]
   end
   return outvals, outvecs, outquns
end

function tsrcalc2(prm,stg,cdo,nf,ctrl,jlist)
   s = ctrl["S"]
   mcalc = ctrl["mcalc"]
   sd = Int(2*s + 1)
   sof = prm[1:15]
   cdf = prmsetter(prm[16:end],stg)
   vtd = Int(ctrl["vtmax"]+1)
   jmin = 0.5*iseven(sd)
   jmax = jlist[end,1]
   jfd = sd*Int(sum(2.0 .* collect(Float64,jmin:jmax) .+ 1.0))
   σcnt = σcount(nf)
   mcd = Int(2*mcalc+(iseven(nf))+1)
   fvls = zeros(Float64,jfd*vtd,σcnt)
   fqns = zeros(Int,jfd*vtd,6,σcnt)
   fvcs = zeros(Float64,Int(sd*(2*jmax+2)*mcd),jfd*vtd,σcnt)
   @time @simd for sc in 1:σcnt
      σ = sc - 1
      mcd = Int(2*mcalc+(σtype(nf,σ)==2)+1)
      tormat, mk, mb = htor(sof,nf,mcalc,σ)
      mcd = Int(2*mcalc+(σtype(nf,σ)==2)+1)
      σt = σtype(nf,σ)
      msd = sd*mcd
      mstrt, mstop = mslimit(nf,mcalc,σ)
      jmsd = Int(mcd*sd*(2*jmax+1))
      jsvd = Int(jfd*vtd)
      jsublist = jlist[isequal.(jlist[:,2],σ), 1] .* 0.5
      for j in jsublist
         jd = Int(2.0*j) + 1
         pull = indpuller(ctrl["vtmax"],mcalc,σt,Int(jd*sd))
         sind, find = jvdest(j,s,ctrl["vtmax"]) 
         tvals, tvecs = tsrdiag(sof,cdf,cdo,tormat,nf,mcalc,mb,mk,j,s,σ,σt)
         fvls[sind:find,sc] = tvals[pull]
         fvcs[1:jd*msd,sind:find,sc] = tvecs[:,pull]
         fqns[sind:find,:,sc] = qngenv(j,s,nf,ctrl["vtmax"],σ)
      end
   end
   return fvls, fvcs, fqns
end
