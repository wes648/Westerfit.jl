

struct μFunc{T <: Number}
   f::FunctionWrapper{SparseMatrixCSC{Float64,Int}, Tuple{RPsi,Int}}
   k::Int
   q::Int
   function μFunc(T::Type,f::Function,k::Int,q::Int)
      new{T}(f,k)
   end
   function μFunc(T::Type,f::Function,k::Int,q::Int)
      new{Float64}(f,k,q)
   end
end

eval_μop(op::OpFunc,ψ::RPsi)::SparseMatrixCSC{Float64,Int} = op.f(ψ,op.k,op.q)
eval_μop(op::OpFunc,ψ::TPsi)::SparseMatrixCSC{Float64,Int} = op.f(ψ,op.k,op.q)
