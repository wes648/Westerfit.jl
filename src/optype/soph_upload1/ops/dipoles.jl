#=
using SparseArrays
Δlist2(J::Real,S::Real)::UnitRange{Int} = Int(abs(J-S)):Int(J+S)
powneg1(k::Real)::Int = isodd(k) ? -1 : 1
struct RPsi
   J::Float64
   S::Float64
   N::UnitRange{Int}
   lng::Int
   function RPsi(J::Number,S::Number)
      J = Float64(J)
      S = Float64(S)
      N = Δlist2(J,S)
      lng = Int((2J+1)*(2S+1))
      new(J,S,N,lng)
   end
   function RPsi(N::Int)
      J = Float64(N)
      N = Δlist2(J,0.0)
      lng = Int(2J+1)
      new(J,0.0,N,lng)
   end
end
function wig3j(j1,j2,j3,m1,m2,m3)::Int
   return iszero(m1+m2+m3)&&(abs(j1-j2)≤j3≤abs(j1+j2))
end
=#

function μ1q_elem(jb,jk,s,nb,kb,nk,kk,q)
   out = powneg1(s+jb-kb+nb+nk)*√((2jb+1)*(2jk+1)*(2nk+1)*(2nb+1))
   out *= wig3j(nb,1,nk, -kb, q, kk)
   out *= wig6j(nk,jk,s, jb,nb, 1)
end

function μkq(ψb::RPsi,ψk::RPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   jb = fill(ψb.J, ψb.lng)
   nb = mapreduce(x->fill(x,2x+1), append!, ψb.N )
   kb = mapreduce(x->-x:x, vcat, ψb.N )
   jk = fill(ψk.J, ψk.lng)'
   nk = mapreduce(x->fill(x,2x+1), append!, ψk.N )'
   kk = mapreduce(x->-x:x, vcat, ψk.N )'
   out = μ1q_elem.(jb,jk,ψk.S,nb,kb,nk,kk, q)
   out = ur(ψb.J,ψb.S) * sparse(out) * ur(ψk.J,ψk.S)
   return dropzeros!(out)
end

μz(ψb::RPsi,ψk::RPsi,k::Int,q::Int) = μkq(ψb,ψk,1,0)
μx(ψb::RPsi,ψk::RPsi,k::Int,q::Int) = √(0.5) .* (-μkq(ψb,ψk,1,1) + μkq(ψb,ψk,1,-1))

function μcos(ψb::TPsi,ψk::TPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   mb = collect(ψb.ms)
   mk = collect(ψk.ms)
   out = abs.(mb' .- mk) .== k
   return sparse(out)
end

μcosα(ψb::TPsi,ψk::TPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int} = μcos(ψb,ψk,k,1)
μcosβ(ψb::TPsi,ψk::TPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int} = μcos(ψb,ψk,k,2)
μcosγ(ψb::TPsi,ψk::TPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int} = μcos(ψb,ψk,k,3)
