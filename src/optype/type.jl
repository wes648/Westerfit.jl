
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
   function TPsi(nf,σ,mc=3)
      if length(nf) > 1
         ms = msgen_indef(nf,mc,σ)
      else
         ms = [msgen(nf,mc,σ)]
         σ = [σ]
         nf = [nf]
      end
      lng = sum(length.(ms))
      new(nf,ms,σ,lng)
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
   f::FunctionWrapper{SparseMatrixCSC{Float64,Int}, Tuple{RPsi,Int}}
   p::Int
   function OpFunc(T::Type,f::Function,p::Int)
      new{T}(f,p)
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

#function tarraysum(a::Array{Int,2},b::Array{Int,2})::Array{Int,2}
#   if size(a,2) < size(b,2)
#      a = hcat(a,zeros(Int,2,size(b,2)-size(a,2)))
#   elseif size(b,2) < size(a,2)
#      b = hcat(b,zeros(Int,2,size(a,2)-size(b,2)))
#   end
#   return a .+ b
#end



#= This section is commented out due to removal of op mutation
#Multiplying an Operator by a number updates the value
function *(v::Number,O::Op)::Op
   out = Op(O.nam,v*O.v,O.rf,O.tp,O.a,O.b,O.b,O.c,O.d)
   return out
end
#This allows a scalar to be distributed among a list of added Operators
function *(v::Number,O::Vector{Op})::Vector{Op}
   out = similar(O)
   for i in 1:length(O)
      out[i] = v*O[i]
   end
   return out
end
#Raising an Operator to a power updates the exponent
function ^(O::Op,n::Int)::Op 
   out = Op(O.nam,O.v,rf=O.rf[:].p.*n,tp=O.tp.*n,
      a=O.a*n,b=O.b*n,c=O.c*n,d=O.d*n)
   return out
end
#function ^(O::Vector{Op},n::Int)::Vector{Op}
#   out = O
#   for i in 2:n
#      out *= O
#   end
#   return out
#end
#Multiplying two Operators multiplies the values & concatenates the powers + functions
function *(O::Op,P::Op)::Op
   Op(O.nam,O.v*P.v,rf=vcat(O.rf, P.rf),tp=tarraysum(O.tp,P.tp),a=O.a+P.a,
      b=O.b+P.b,c=O.c+P.c,d=O.d+P.d)
end
function *(O::Op,P::Vector{Op})::Vector{Op}
   OP = similar(P)
   for i in eachindex(P)
      OP[i] = O * P[i]
   end
   return OP
end
function *(O::Vector{Op},P::Op)::Vector{Op}
   OP = similar(O)
   for i in eachindex(O)
      OP[i] = O[i] * P
   end
   return OP
end
function *(O::Vector{Op},P::Vector{Op})::Vector{Op}
   OP = Vector{Op}(undef,length(O)+length(P))
   for i in eachindex(OP)
      o = floor(Int,(i-1)/length(O)) + 1
      p = mod(i-1,length(P)) + 1
      OP[i] = O[o] * P[p]
   end
   return OP
end
#Adding two Operators generates a vector of Operators
+(O::Op,P::Op)::Vector{Op} = vcat(O,P)
+(O::Op,P::Vector{Op})::Vector{Op} = vcat(O,P)
+(O::Vector{Op},P::Op)::Vector{Op} = vcat(O,P)
+(O::Vector{Op},P::Vector{Op})::Vector{Op} = vcat(O,P)
#Subtraction rules
-(O::Op,P::Op)::Vector{Op} = vcat(O,-1*P)
-(O::Op,P::Vector{Op})::Vector{Op} = vcat(O,-1*P)
-(O::Vector{Op},P::Op)::Vector{Op} = vcat(O,-1*P)
-(O::Vector{Op},P::Vector{Op})::Vector{Op} = vcat(O,-1*P)
#Δlist(J::Real,S::Real)::Array{Int} = collect(Int(abs(J-S)):Int(J+S))
#function ngen(J::Real,S::Real)::Vector{Int}
#   ns = Δlist(J,S)
#   out = fill(ns[1],2*ns[1]+1)
#   for i in 2:length(ns)
#      out = vcat(out,fill(ns[i],2*ns[i]+1))
#   end
#   return out
#end
#function kgen(J::Real,S::Real)::Vector{Int}
#   ns = Δlist(J,S)
#   out = collect(-ns[1]:ns[1])
#   for i in 2:length(ns)
#      n = ns[i]
#      out = vcat(out,collect(-n:n))
#   end
#   return out
#end

struct Psi
   J::Float64
   S::Float64
   N::UnitRange{Int}
   K::Vector{UnitRange{Int}}
   nf::Vector{Int}
   ms::Vector{StepRange{Int,Int}}
   σ::Vector{Int}
   lng::Int
   #Psi(J::Real,S::Real) = new(Float64(J),Float64(S),ngen(J,S),kgen(J,S),Int((2J+1)*(2S+1)))
   #Psi(n::Int) = Psi(Float64(n),0.0)
   function Psi(J::Number=0,S::Number=0;nf=0,σ=0,mc=0) #spin (torsion) rotation
      J = convert(Float64,J)
      S = convert(Float64,S)
      N = Δlist2(J,S)
      K = kgen(N)
      if length(nf) > 1
         ms = msgen_indef(nf,mc,σ)
      else
         ms = [msgen(nf,mc,σ)]
         nf = [nf]
         σ = [σ]
      end
      lng = convert(Int,(2J+1)*(2S+1))
      new(J,S,N,K,nf,ms,σ,lng)
   end
   function Psi(N::Int=0;nf=0,σ=0,mc=3) # (torsion) rotation
      J = convert(Float64,N)
      S = 0.0
      N = ngen(J,S)
      K = kgen(J,S)
      if length(nf) > 1
         ms = msgen_indef(nf,mc,σ)
      else
         ms = msgen(nf,mc,σ)
         nf = [nf]
         σ = [σ]
      end
      lng = convert(Int,(2J+1)*(2S+1))
      new(J,S,N,K,nf,ms,σ,lng)
   end
   function Psi(nf=0,σ=0,mc=3) #torsion
      J = 0.0
      S = 0.0
      N = Δlist2(0,0)
      K = kgen(0:0)
      if length(nf) > 1
         ms = msgen_indef(nf,mc,σ)
      else
         ms = [msgen(nf,mc,σ)]
         σ = [σ]
         nf = [nf]
      end
      lng = convert(Int,(2J+1)*(2S+1))
      new(J,S,N,K,nf,ms,σ,lng)
   end
end
#mutable struct Op_old
#   #v for value as in parameter value
#   v::Float64
#   #rp for rotation Operator powers
#   rp::Vector{Int}
#   #rf for rotation Operator functions
#   rf::Vector{Function}
#   #tor powers. each column is for the consecutive rotors 
#   tp::Array{Int,2}
#   #a,b,c,d are powers for N^2a (NS)^b S^2c Nz^d
#   #this is to initialize the purely diagonal part of the matrix first
#   a::Int
#   b::Int
#   c::Int
#   d::Int
#   #forces the vector structure for even a single function
#   #Op(v::Number,p::Int,f::Function;a=0,b=0,c=0,d=0) = new(Float64(v),[p],[f],a,b,c,d)
#   function Op_old(v::Number;rp=Vector{Int}[],rf=Vector{Function}[],
#               tp=zeros(Int,2,0),a=0,b=0,c=0,d=0)
#      return new(Float64(v),rp,rf,tp,a,b,c,d)
#   end
#   Op_old(O::Op) = new(Float64(O.v),O.rp,O.rf,O.tp,O.a,O.b,O.c,O.d)
#end
=#
