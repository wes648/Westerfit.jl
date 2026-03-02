
Ψ = RPsi(1)

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

nzob = Op("BK",[OpFunc(Nz,2)])
H_test = [Op("BK",[OpFunc(Nz,2)]);
          Op("BN",[OpFunc(N2,1)]);
          Op("B±",[OpFunc(Npm,2)])]



#= 
H = spzeros(3,3)
for i ∈ 1:3
	H += eval_rop(H_test[i].rf[1], Ψ)
end

=#