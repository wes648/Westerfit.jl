
using LinearAlgebra
using SparseArrays
import Base.*
import Base.^
import Base.+
import Base.-

struct psi
   N::Int
end

struct op
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
*(v::Number,O::op)::op = op(O.v*v,O.p,O.f)
#Raising an operator to a power updates the exponent
^(O::op,n::Int)::op = op(O.v,O.p*n,O.f)
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
      OP[i] = O * P[i]
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
function *(O::op,ψ::psi)
   out = O.v
   for i in eachindex(O.p)
      out *= (O.f[i](ψ))^O.p[i]
   end
   return out
end
#This allows the basis set to be distributed among a list of added operators
function *(O::Vector{op},ψ::psi)
   out = O[1]*ψ
   for i in 2:length(O)
      out += O[i]*ψ
   end
   return out
end
#This allows a scalar to be distributed among a list of added operators
function *(v::Number,O::Vector{op})::Vector{op}
   out = v*O[1]
   for i in 2:length(O)
      out += v*O[i]
   end
   return out
end

#Get it? \squarert? It's a more stable variant of sqrt that defaults to zero when
#   given a negative input value
function □rt(x)::Float64
   if x ≤ zero(x)
      return 0.0
   else
      return √x
   end
end
function fh(x::Int,y::Int)::Float64
   out = □rt((x-y)*(x+y+1))
   return out
end

#The four nice angular momentum operators
function nz(ψ::psi)::Diagonal{Float64, Vector{Float64}}
   Diagonal(collect(Float64,-ψ.N:ψ.N))
end
function nt(ψ::psi)::Diagonal{Float64, Vector{Float64}}
   Diagonal(fill(√(ψ.N*(ψ.N+1)),2*ψ.N +1))
end
function np(ψ::psi)::SparseMatrixCSC{Float64, Int64}
   out = fh.(ψ.N,collect(-ψ.N:ψ.N-1))
   out = spdiagm(-1=>out)
   return out
end
function nm(ψ::psi)::SparseMatrixCSC{Float64, Int64}
   out = fh.(ψ.N,collect(-ψ.N:ψ.N-1))
   out = spdiagm(1=>out)
   return out
end

### TEST BIT FOR DEMO ###
using Symbolics
@variables A B C
A = 3000.
B = 1500.
C = 1000.
Nz = op(1.0,1,nz)
Nt = op(1.0,1,nt)
Np = op(1.0,1,np)
Nm = op(1.0,1,nm)
Nx = 0.5*(Np + Nm)
Ny = 0.5*im*(Np-Nm)
#Setting the value of psi defines N
ϕ = psi(1)

#The effective vs pretty forms of the asymmetric top hamiltonian
H = (A - 0.5*(B+C))*Nz^2 + 0.5*(B+C)*Nt^2 + 0.25*(B-C)*(Np^2 + Nm^2)
Hv2 = A*Nz^2 + B*Nx^2 + C*Ny^2

H*ϕ
Hv2*ϕ


