

struct μFunc
   μv::Float64
   fr::FunctionWrapper{SparseMatrixCSC{Float64,Int}, Tuple{RPsi,Int}}
   ft::FunctionWrapper{SparseMatrixCSC{Float64,Int}, Tuple{TPsi,Int}}
   k::Int
   q::Int
   function μFunc(μv::Float64,fr::Function,ft::Function,k::Int,q::Int)
      new(μv,fr,ft,k,q)
   end
end

function eval_μop_r(op::μFunc,ψb::RPsi,ψk::RPsi)::SparseMatrixCSC{Float64,Int} 
   op.fr(ψb,ψk,op.k,op.q)
end
function eval_μop_t(op::μFunc,ψb::TPsi,ψk::TPsi)::SparseMatrixCSC{Float64,Int} 
   op.ft(ψb,ψk,op.k,op.q)
end

function int_enact(μp,μf,ψb,ψk,bvec,kvec)::SparseMatrixCSC{Float64,Int}
   out = kron(μp,eval_μop_t(μf,ψb.T,ψk.T),eval_μop_r(μf,ψb.R,ψk.R))
end

function jbjklister(jmin,jmax,mΔj)::Matrix{Float64}
   out = zeros(0,2)
   for j in jmin:jmax
      p1 = collect(max(jmin,(j-mΔj)):(j+mΔj))
      out = vcat(out,[fill(j,length(p1)) p1])
   end
   return out
end



################################################################################
###     Parametric Typing Testing
#=
using FunctionWrappers
import FunctionWrappers: FunctionWrapper

struct StabOp{T<:Number}
   f::FunctionWrapper{T, Tuple{Float64,Int}}
   p::Int
   function StabOp(T::Type,f::Function,p::Int)
      new{T}(f,p)
   end
end

function eval_op(op::StabOp{T},arg)::T where {T<:Number}
   op.f(arg, op.p)
end

e(x::Float64,p::Int) = x^p
f(x::Float64,p::Int) = im * (x^p)
g(x::Float64,p::Int) = cos(p*x)
h(x::Float64,p::Int) = im*sin(p*x)

E = StabOp(Float64,e,2)
F = StabOp(ComplexF64,f,3)
G = StabOp(Float64,g,4)
H = StabOp(ComplexF64,h,5)
=#