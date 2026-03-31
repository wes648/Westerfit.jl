
eh(x::Real)::Float64 = x*(x+1)
function n2(ns::UnitRange{Int},p::Int)::Vector{Float64}
   mapreduce(x->fill(eh(x)^p, 2x+1), vcat, ns)
end
function N2(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int64}
	return spdiagm(n2(ψ.N,p))
end

nz(ks::Vector{UnitRange{Int}},p::Int)::Vector{Float64} = map(x->x^p, reduce(vcat, ks))
function Nz(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int64}
	return spdiagm(nz(ψ.K,p))
end

□rt(x::Real)::Float64 =√(x*(x>zero(x)))
fh(x::Real,y::Real)::Float64 = □rt((x-y)*(x+y+1))
function np(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int64}
   if p ≠ 0 &&p ≤ ψ.lng
      ns = reduce(vcat, [fill(n,2n+1) for n ∈ ψ.N])[1+p:end]
      ks = reduce(vcat, ψ.K)[1+p:end]
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
	return tplus!(np(ψ,p))
end


function p_tor(ψ::TPsi,p::Int,tid::Int)::SparseMatrixCSC{Float64, Int}
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
   #@show tid
   #torsetter!(ψ,tid,out)
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

Pα(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = p_tor(ψ,p,1)
Pβ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = p_tor(ψ,p,2)
Pγ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = p_tor(ψ,p,3)

cosα(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = cos_tor(ψ,p,1)
cosβ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = cos_tor(ψ,p,2)
cosγ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = cos_tor(ψ,p,3)
vncα(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = vnc_tor(ψ,p,1)
vncβ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = vnc_tor(ψ,p,2)
vncγ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = vnc_tor(ψ,p,3)

sinα(ψ::TTPsi,p::Int,q::Int)::SparseMatrixCSC{ComplexF64, Int} = sin_tor(ψ,p,1)
sinβ(ψ::TTPsi,p::Int,q::Int)::SparseMatrixCSC{ComplexF64, Int} = sin_tor(ψ,p,2)
sinγ(ψ::TTPsi,p::Int,q::Int)::SparseMatrixCSC{ComplexF64, Int} = sin_tor(ψ,p,3)

