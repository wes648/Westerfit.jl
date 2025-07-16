#c1 testing
using FunctionWrappers
import FunctionWrappers: FunctionWrapper
using LinearAlgebra
using BenchmarkTools

return_type(f::Function)::Type = eltype(Base.return_types(f)[1])

function freal(l::Int,p::Int)::Array{Float64,2}
   return rand(Float64,l,l)^p
end
function gcomp(l::Int,p::Int)::Array{ComplexF64,2}
   return rand(ComplexF64,l,l)^p
end

struct OpFunc{T <: Number}
   f::FunctionWrapper{Array{T,2}, Tuple{Int,Int}}
   p::Int
   function OpFunc(T::Type,f::Function,p::Int)
      new{T}(f,p)
   end
   function OpFunc(f::Function,p::Int)
      new{return_type(f)}(f,p)
   end
end


struct OpFunc{T <: Number}
   f::FunctionWrapper{Array{T,2}, Tuple{Psi,Int}}
   p::Int
   function OpFunc(T::Type,f::Function,p::Int)
      new{T}(f,p)
   end
   function OpFunc(f::Function,p::Int)
      new{return_type(f)}(f,p)
   end
end


function enact(l::Int,f::OpFunc{T})::Array{T,2} where {T<:Number}
   return f.f(l,f.p)
end
function enact2(init::Array{T,2},l::Int,f::OpFunc)::Array{T,2} where {T<:Number}
   return init + f.f(l,f.p)
end


fop = OpFunc(freal,2)
gop = OpFunc(gcomp,2)

