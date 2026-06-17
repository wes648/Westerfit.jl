
using LinearAlgebra, WIGXJPFjl, SparseArrays

powneg1(k::Real)::Int = isodd(k) ? -1 : 1

function adiagin(a::AbstractMatrix)::StepRange{Int,Int}
   @assert size(a,1)==size(a,2)
   l = size(a,1)
   return l:(l-1):(l*(l-1)+1)
end

function ur(n::Int)::SparseMatrixCSC{Float64, Int}
   if !iszero(n)
   out = spzeros(1:2n+1,1:2n+1)
   out[diagind(out)[1:n]] .= -√.5
   out[diagind(out)[n+2:end]] .= √.5
   out[adiagin(out)] .= √.5
   out[n+1,n+1] = 1.0
   return sparse(out)
   else
      return sparse(1.0I,1,1)
   end
end

function ur2(n::Int)::SparseMatrixCSC{Float64, Int}
   if !iszero(n)
   out = spzeros(1:2n+1,1:2n+1)
   out[diagind(out)[1:n]] .= √.5 .*powneg1.(1:n)
   out[diagind(out)[n+2:end]] .= √.5 
   out[adiagin(out)[1:n]] .= √.5
   out[adiagin(out)[n+2:end]] .= -√.5 .*powneg1.(1:n)
   out[n+1,n+1] = 1.0
   return sparse(out)
   else
      return sparse(1.0I,1,1)
   end
end


function T(l::Int,q::Int,n::Int)::SparseMatrixCSC{Float64,Int}
   ns = fill(n,2n+1)
   ks = collect(-n:n)
   out = @. wig3j(ns',l,ns, -ks',q,ks) * (-1)^(ns' - ks')
   return out
end
function R(l::Int,q::Int,n::Int)::SparseMatrixCSC{Float64,Int}
   ns = fill(n,2n+1)
   ks = collect(-n:n)
   out = @. wig3j(ns',l,ns, -ks',q,ks) * (-1)^(ns' - ks')
   out = ur(n)*sparse(out)*ur(n)
   return dropzeros(out)
end

function H(l::Int,q::Int,n::Int)::SparseMatrixCSC{Float64,Int}
   out = T(l,abs(q),n) + (-1)^q .* T(l,-abs(q),n)
   return dropzeros(out)
end

function W(l::Int,q::Int,n::Int)::SparseMatrixCSC{Float64,Int}
   out = T(l,abs(q),n) + (-1)^q .* T(l,-abs(q),n)
   out = ur(n)*out*ur(n)
   return dropzeros(out)
end