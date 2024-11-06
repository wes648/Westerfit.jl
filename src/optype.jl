
using LinearAlgebra
using SparseArrays
using BenchmarkTools
import Base.*
import Base.^
import Base.+
import Base.-


tplus!(a::Diagonal)::Diagonal = a
tplus!(a::Array{Float64,2})::Array{Float64,2} = hermitianpart!(2a)
function tplus!(a::SparseMatrixCSC{Float64, Int64})::SparseMatrixCSC{Float64, Int64}
   a .+= permutedims(a)
end

struct psi
   N::Int
   K::Vector{Int}
   lng::Int
   psi(N::Int) = new(N,collect(-N:N),2N+1)
end

mutable struct op
   #v for value as in parameter value
   v::Number
   #p for power as in operator power
   p::Vector{Int64}
   #f for function as in the matrix element function
   f::Vector{Function}
   #forces the vector structure for even a single function
   op(v::Number,p::Int,f::Function) = new(v,[p],[f])
   op(v::Number,p::Vector{Int64},f::Vector{Function}) = new(v,p,f)
end

#Multiplying an operator by a number updates the value
*(v::Number,O::op)::op = op(v*O.v,O.p,O.f)
#Raising an operator to a power updates the exponent
^(O::op,n::Int)::op = op(O.v,n*O.p,O.f)

function ^(O::Vector{op},n::Int)::Vector{op}
   out = O
   for i in 2:n
      out *= O
   end
   return out
end
#Multiplying two operators multiplies the values & concatenates the powers + functions
*(O::op,P::op)::op = op(O.v*P.v,vcat(O.p, P.p),vcat(O.f, P.f))
function *(O::op,P::Vector{op})::Vector{op}
   OP = similar(P)
   for i in eachindex(P)
      OP[i] = O * P[i]
   end
   return OP
end
function *(O::Vector{op},P::op)::Vector{op}
   OP = similar(O)
   for i in eachindex(O)
      OP[i] = O[i] * P
   end
   return OP
end
function *(O::Vector{op},P::Vector{op})::Vector{op}
   OP = Vector{op}(undef,length(O)+length(P))
   for i in eachindex(OP)
      o = floor(Int,(i-1)/length(O)) + 1
      p = mod(i-1,length(P)) + 1
      OP[i] = O[o] * P[p]
   end
   return OP
end

#Adding two operators generates a vector of operators
+(O::op,P::op)::Vector{op} = vcat(O,P)
+(O::op,P::Vector{op})::Vector{op} = vcat(O,P)
+(O::Vector{op},P::op)::Vector{op} = vcat(O,P)
+(O::Vector{op},P::Vector{op})::Vector{op} = vcat(O,P)
#Subtraction rules
-(O::op,P::op)::Vector{op} = vcat(O,-1*P)
-(O::op,P::Vector{op})::Vector{op} = vcat(O,-1*P)
-(O::Vector{op},P::op)::Vector{op} = vcat(O,-1*P)
-(O::Vector{op},P::Vector{op})::Vector{op} = vcat(O,-1*P)


function ur(n::Int)::SparseMatrixCSC{Float64, Int64}
   out = Diagonal(vcat(fill(-√.5,n), 1.0, fill(√.5,n)))
   out += rotl90(Diagonal(vcat(fill(√.5,n), 0.0, fill(√.5,n))))
   return sparse(out)
end
#*(O::op,ψ::psi) = O.v * (O.f(ψ))^O.p
#Multiplying an operator by the wavefunction gnerates the matrix
#   Yeah this isn't the most QMly sound notation but it should be
#   fairly comfortable & easy for spectroscopists to follow
function *(O::op,ψ::psi)#::SparseMatrixCSC{Float64, Int64}
   out = O.v*I(ψ.lng)
   for i in eachindex(O.p)
      out *= (O.f[i](ψ,O.p[i]))
   end
   tplus!(out)
   return out
end
#This allows the basis set to be distributed among a list of added operators
function *(O::Vector{op},ψ::psi)::SparseMatrixCSC{Float64, Int64}
   out = O[1]*ψ
   for i in 2:length(O)
      out += O[i]*ψ
   end
   return out
end
#This allows a scalar to be distributed among a list of added operators
function *(v::Number,O::Vector{op})::Vector{op}
   out = similar(O)
   for i in 1:length(O)
      out[i] = v*O[i]
   end
   return out
end

#Get it? \squarert? It's a more stable variant of sqrt that defaults to zero when
#   given a negative input value
eh(x::Real)::Float64 = √(x*(x+1))
□rt(x::Real)::Float64 =√(x*(x>zero(x)))
fh(x::Real,y::Real)::Float64 = □rt((x-y)*(x+y+1))

#The four nice angular momentum operators
function nz(ψ::psi,p::Int)::Diagonal{Float64, Vector{Float64}}
   Diagonal(ψ.K)^p
end
function nt(ψ::psi,p::Int)::Diagonal{Float64, Vector{Float64}}
   Diagonal(fill(eh(ψ.N)^p, ψ.lng))
end
function np(ψ::psi,p::Int)::SparseMatrixCSC{Float64, Int64}
   ns = fill(ψ.N,ψ.lng-p)
   part = ones(length(ns))
   if p ≤ length(ns) && p ≠ 0
      ks = ψ.K[1+p:end]
      part = ones(length(ks))
      for o in 1:p
         part .*= fh.(ns,ks.-o)
      end
   #end#original if
      out = spzeros(ψ.lng,ψ.lng)
      out[diagind(out,-p)] = part
   elseif p > length(ns) && p ≠ 0
      out = spzeros(ψ.lng,ψ.lng)
      out[diagind(out,-p)] = part
   else
      out = spdiagm(ones(ψ.lng))
   end
   return out
end


### TEST BIT FOR DEMO ###
Nz = op(1.0,1,nz)
Nt = op(1.0,1,nt)
Np = op(1.0,1,np)
#Nm = op(1.0,1,nm)
Nz2 = op(1.0,2,nz)
Nt2 = op(1.0,2,nt)
Np2 = op(1.0,2,np)
#Nm2 = op(1.0,2,nm)
ψ = psi(3)
#Setting the value of psi defines N
#ϕ = psi(1)

#The effective vs pretty forms of the asymmetric top hamiltonian
#H = (A - 0.5*(B+C))*Nz^2 + 0.5*(B+C)*Nt^2 + 0.25*(B-C)*(Np^2 + Nm^2)
BK = 1.75
BN = 1.25
Bp = 0.125
H = BK*Nz2 + BN*Nt2 + Bp*(Np2)
#Hv2 = A*Nz^2 + B*Nx^2 + C*Ny^2

#H*ϕ
#Hv2*ϕ


