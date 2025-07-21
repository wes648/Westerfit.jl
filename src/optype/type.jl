
import Base.*
import Base.^
import Base.+
import Base.-
import Base: convert

struct RPsi
   J::Float64
   S::Float64
   N::UnitRange{Int}
   K::Vector{UnitRange{Int}}
   lng::Int
   function RPsi(J::Number,S::Number)
      J = convert(Float64,J)
      S = convert(Float64,S)
      N = Δlist2(J,S)
      K = kgen(N)
      lng = convert(Int,(2J+1)*(2S+1))
      new(J,S,N,K,lng)
   end
   function RPsi(N::Int)
      J = convert(Float64,N)
      S = 0.0
      N = Δlist2(J,0.0)
      K = kgen(N)
      lng = convert(Int,(2J+1)*(2S+1))
      new(J,S,N,K,lng)
   end
end
struct TPsi
   nf::Vector{Int}
   ms::Vector{StepRange{Int,Int}}
   σ::Vector{Int}
   lng::Int
   function TPsi(nf::Vector{Int},σ::Vector{Int},mc=3)
      if length(nf) > 1
         ms = msgen(nf,mc,σ)
      else
         ms = [msgen(nf,mc,σ[1])]
      #   σ = [σ]
         nf = nf
      end
      lng = sum(length.(ms))
      new(nf, ms, σ, lng)
   end
   function TPsi(nf::Int,σ::Int,mc=3)
      new([nf],[σ],mc,2mc+1)
   end
end
struct VPsi
   #each column refer to a vibrational polyad state
   #each row refers to the quanta of a given mode
   #thus nus[a,b] is the quanta of mode a in state b
   nus::Matrix{Int}
end
mutable struct Psi
   R::RPsi
   T::TPsi
end

convert(T::Type{Psi},ϕ::RPsi) = Psi(ϕ,TPsi(0,0,0))
convert(T::Type{Psi},ϕ::TPsi) = Psi(RPsi(0),ϕ)

struct OpFunc{T <: Number}
   f::FunctionWrapper{SparseMatrixCSC{T,Int}, Tuple{RPsi,Int}}
   p::Int
   function OpFunc(T::Type,f::Function,p::Int)
      new{T}(f,p)
   end
   function OpFunc(f::Function,p::Int)
      new{Float64}(f,p)
   end
end
#struct OpFunc
#   f::FunctionWrapper{SparseMatrixCSC{Float64,Int}, Tuple{RPsi,Int}}
#   p::Int
#end
eval_rop(op::OpFunc,ψ::RPsi)::SparseMatrixCSC{Float64,Int} = op.f(ψ,op.p)
struct Op
   #this is a nam field
   nam::String
   #v for value as in parameter value
   #v::Float64
   #rf for rotation Operators
   rf::Vector{OpFunc}
   #tor powers. each column is for the consecutive rotors 
   tp::Array{Int,2}
   #a,b,c,d are powers for N^2a (NS)^b S^2c Nz^d
   #this is to initialize the purely diagonal part of the matrix first
   a::Int
   b::Int
   c::Int
   d::Int
   #forces the vector structure for even a single function
   function Op(nam::String;rf=Vector{OpFunc}[],tp=zeros(Int,2,0),a=0,b=0,c=0,d=0)
      return new(nam,rf,tp,a,b,c,d)
   end
   function Op(nam::String,rf::Vector{OpFunc},abcd::Vector{Int},tp::Array{Int,2})
      return new(nam,rf,tp,abcd[1],abcd[2],abcd[3],abcd[4])
   end
   Op(O::Op) = new(O.nam,O.rf,O.tp,O.a,O.b,O.c,O.d)
end

