
eh(x::Real)::Float64 = x*(x+1)
□rt(x::Real)::Float64 =√(x*(x>zero(x)))
fh(x::Real,y::Real)::Float64 = □rt((x-y)*(x+y+1))
fhv(x::Float64,y::Int)::Float64 = □rt(x - eh(y))

function nz(ns::UnitRange{Int},p::Int)::Vector{Float64}
   mapreduce(x->(-x:x).^p, append!, ns)
end

function n2(ns::UnitRange{Int},p::Int)::Vector{Float64}
   mapreduce(x->fill(eh(x)^p, 2x+1), append!, ns)
end

function N2(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int64}
	return spdiagm(n2(ψ.N,p))
end

function Nz(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int64}
	return spdiagm(nz(ψ.N,p))
end

function Np(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int64}
   if p ≠ 0 &&p ≤ ψ.lng
      ns = reduce(vcat, [fill(n,2n+1) for n ∈ ψ.N])[1+p:end]
      ks = mapreduce(x->-x:x, vcat, ψ.N)[1+p:end]
      out = ones(length(ks))
      out = prod(fh.(ns,ks .- collect(1:p)'),dims=2)[:]
      out = spdiagm(-p=>out)
   elseif p ≠ 0 && p > ψ.lng
      out = spzeros(ψ.lng,ψ.lng)
   else
      out = spdiagm(ones(ψ.lng))
   end
   return out
end
function Npm(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int64}
	return tplus!(Np(ψ,p))
end


function s2(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int}
   return spdiagm( fill(eh(ψ.S)^p,ψ.lng) )
end

function skqpart(j,s,k,q,nb,kb,nk,kk)::Float64
   return wig3j(nb,k,nk,-kb,q,kk)*wig6j(s,nb,j,nk,s,k)*jnred(nb,nk)*powneg1(-kb)
end
function s0part(j,s,nb,kb,nk,kk)::Float64
   return skqpart(j,s,1,0,nb,kb,nk,kk)::Float64
end

function Sz(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   #pr = [T0_0 T2_0 T2_1 T2_2]
   J = ψ.J
   S = ψ.S
   nds = nindsgen(ψ.N)
   out = spzeros(ψ.lng,ψ.lng)
   sfact = nred(S)*powneg1(J+S)
   for i ∈ 1:length(ψ.N), j ∈ max(i-1,1):min(i+1,length(ψ.N))
      nb = ψ.N[j]; nk = ψ.N[i]; Δ = nk - nb
      blck = view(out,nds[j],nds[i])
      frac = wig6j(S,nb,J,nk,S,1)*jnred(nb,nk)*sfact
      #p = Δ
      dest = diagind(blck,Δ)
      kl = (-nk:nk)[(1:length(dest)).+δi(1,Δ)]
      blck[dest] = srelem(frac,nb,nk,kl,1,0)
   end
   return dropzeros!(out^p)
end
#sq(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int} = sq_op(ψ.J,ψ.S,ψ.N,ψ.K,ψ.lng,p)

function Sq(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   J = ψ.J
   S = ψ.S
   nds = nindsgen(ψ.N)
   out = spzeros(ψ.lng,ψ.lng)
   sfact = nred(S)*powneg1(J+S+p)
   for i ∈ 1:length(ψ.N), j ∈ max(i-p,1):min(i+p,length(ψ.N))
      nb = ψ.N[j]; nk = ψ.N[i]; Δ = nk - nb
      blck = view(out,nds[j],nds[i])
      frac = wig6j(S,nb,J,nk,S,p)*jnred(nb,nk)*sfact
      p = Δ+p
      dest = diagind(blck,p)
      kl = (-nk:nk)[(1:length(dest))] .+ Δ
      blck[dest] = srelem(frac,nb,nk,kl,p,-p)
   end
   dropzeros!(out)
   return out
end


Sp(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int} = p*Sq(ψ,p,0)
Sm(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int} = p*Sq(ψ,-p,0)
Spm(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int} = tplus!(p*Sq(ψ,p,0))
function Sx(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int}
#   out  = sq_op(ψ.J,ψ.S,ψ.N,ψ.K,ψ.lng, 1)
#   out -= sq_op(ψ.J,ψ.S,ψ.N,ψ.K,ψ.lng,-1)
   out  = Sq(ψ, 1,0)
   out -= Sq(ψ,-1,0)
   out *= -√.5
   out ^= p
   return out
end


function Pt(ψ::TPsi,p::Int,tid::Int)::SparseMatrixCSC{Float64, Int}
   if iszero(ψ.σ)
      out = map(x -> abs(x)^p, ψ.ms)
      if iseven(p)
         out = sparse(1:ψ.l, 1:ψ.l, out)
      else isodd(p)
         out = sparse(1:ψ.l, ψ.l:-1:1, out)
      end
   else
      out = sparse(1:ψ.l, 1:ψ.l, ψ.ms .^p)
   end
   return out
end

function cos_tor(ψ::TPsi,p::Int,tid::Int)::SparseMatrixCSC{Float64, Int}
   p = floor(Int, p/(ψ.nf * (1+iseven(ψ.nf)) ))
   out = spdiagm(p=>fill(0.5,ψ.l-p),-p=>fill(0.5,ψ.l-p))
   if iszero(ψ.σ)
      u = ul(ψ.l)
      out = dropzeros!(sand(out,u))
   end
   #torsetter!(ψ,tid,out)
   return out
end
function vnc_tor(ψ::TPsi,p::Int,tid::Int)::SparseMatrixCSC{Float64, Int}
   p = floor(Int, p/(ψ.nf * (1+iseven(ψ.nf)) ))
   out = spdiagm(0=>ones(ψ.l),p=>fill(-0.5,ψ.l-p),-p=>fill(-0.5,ψ.l-p))
   if iszero(ψ.σ)
      u = ul(ψ.l)
      out = dropzeros!(sand(out,u))
   end
   #torsetter!(ψ,tid,out)
   return out
end

function sin_tor(ψ::TTPsi,p::Int,tid::Int)::SparseMatrixCSC{ComplexF64, Int}
   l = 2ψ.mc+1
   out = spdiagm(p=>fill(0.5im,l-p),-p=>fill(-0.5im,l-p))
   if iszero(ψ.σs[tid])
      u = ur(ψ.mc)
      out = dropzeros!(u*out*u)
   end
   torsetter!(ψ,tid,out)
   return out
end

Pα(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = Pt(ψ,p,1)
Pβ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = Pt(ψ,p,2)
Pγ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = Pt(ψ,p,3)

cosα(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = cos_tor(ψ,p,1)
cosβ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = cos_tor(ψ,p,2)
cosγ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = cos_tor(ψ,p,3)
vncα(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = vnc_tor(ψ,p,1)
vncβ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = vnc_tor(ψ,p,2)
vncγ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = vnc_tor(ψ,p,3)

sinα(ψ::TTPsi,p::Int,q::Int)::SparseMatrixCSC{ComplexF64, Int} = sin_tor(ψ,p,1)
sinβ(ψ::TTPsi,p::Int,q::Int)::SparseMatrixCSC{ComplexF64, Int} = sin_tor(ψ,p,2)
sinγ(ψ::TTPsi,p::Int,q::Int)::SparseMatrixCSC{ComplexF64, Int} = sin_tor(ψ,p,3)

