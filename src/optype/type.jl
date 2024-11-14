
using JET
using LinearAlgebra
using SparseArrays
using BenchmarkTools
import Base.*
import Base.^
import Base.+
import Base.-

Δlist(J::Real,S::Real)::Array{Int} = collect(Int(abs(J-S)):Int(J+S))
function ngen(J::Real,S::Real)::Array{Int}
   ns = Δlist(J,S)
   out = fill(ns[1],2*ns[1]+1)
   for i in 2:length(ns)
      out = vcat(out,fill(ns[i],2*ns[i]+1))
   end
   return out
end
function kgen(J::Real,S::Real)::Array{Int}
   ns = Δlist(J,S)
   out = collect(-ns[1]:ns[1])
   for i in 2:length(ns)
      n = ns[i]
      out = vcat(out,collect(-n:n))
   end
   return out
end
#tplus!(a::Diagonal)::SparseMatrixCSC{Float64, Int} = sparse(a)
#tplus!(a::Array{Float64,2})::Array{Float64,2} = hermitianpart!(2a)
function tplus!(a::SparseMatrixCSC{Float64, Int})::SparseMatrixCSC{Float64, Int}
   a .+= permutedims(a)
end

struct psi
   S::Float64
   J::Float64
   N::Vector{Int}
   K::Vector{Int}
   nf::Vector{Int}
   ms::Array{Int,2}
   lng::Int
   psi(N::Int) = new(Float64(N),0.0,fill(N,2N+1),collect(-N:N),2N+1)
   psi(J::Real,S::Real) = new(Float64(J),Float64(S),ngen(J,S),kgen(J,S),Int((2J+1)*(2S+1)))
end

mutable struct op_old
   #v for value as in parameter value
   v::Float64
   #p for power as in operator power
   rp::Vector{Int}
   #f for function as in the matrix element function
   f::Vector{Function}
   #forces the vector structure for even a single function
   op_old(v::Number,rp::Int,f::Function) = new(Float64(v),[rp],[f])
   op_old(v::Number,rp::Vector{Int},f::Vector{Function}) = new(Float64(v),rp,f)
end
mutable struct op
   #v for value as in parameter value
   v::Float64
   #rp for rotation operator powers
   rp::Vector{Int}
   #rf for rotation operator functions
   rf::Vector{Function}
   #tor powers. each column is for the consecutive rotors 
   tp::Array{Int,2}
   #a,b,c,d are powers for N^2a (NS)^b S^2c Nz^d
   #this is to initialize the purely diagonal part of the matrix first
   a::Int
   b::Int
   c::Int
   d::Int
   #forces the vector structure for even a single function
   #op(v::Number,p::Int,f::Function;a=0,b=0,c=0,d=0) = new(Float64(v),[p],[f],a,b,c,d)
   function op(v::Number;rp=Vector{Int}[],rf=Vector{Function}[],
               tp=zeros(Int,0,2),a=0,b=0,c=0,d=0)
      return new(Float64(v),rp,rf,tp,a,b,c,d)
   end
   op(O::op) = new(Float64(O.v),O.rp,O.rf,O.tp,O.a,O.b,O.c,O.d)
end

#Multiplying an operator by a number updates the value
*(v::Number,O::op)::op = op(v*O.v,O.p,O.f)
#Raising an operator to a power updates the exponent
function ^(O::op,n::Int)::op 
   out = op(O)
   out.rp .*= n
   out.tp .*= n
   out.a *= n
   out.b *= n
   out.c *= n
   out.d *= n
   return out
end

function ^(O::Vector{op},n::Int)::Vector{op}
   out = O
   for i in 2:n
      out *= O
   end
   return out
end
#Multiplying two operators multiplies the values & concatenates the powers + functions
function *(O::op,P::op)::op
   op(O.v*P.v,rp=vcat(O.rp,P.rp),rf=vcat(O.rf, P.rf),tp=(O.tp .+ P.tp),a=O.a+P.a,
      b=O.b+P.b,c=O.c+P.c,d=O.d+P.d)
end
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


function ur(n::Int)::SparseMatrixCSC{Float64, Int}
   out = Diagonal(vcat(fill(-√.5,n), 1.0, fill(√.5,n)))
   out += rotl90(Diagonal(vcat(fill(√.5,n), 0.0, fill(√.5,n))))
   return sparse(out)
end
#*(O::op,ψ::psi) = O.v * (O.f(ψ))^O.p
#Multiplying an operator by the wavefunction gnerates the matrix
#   Yeah this isn't the most QMly sound notation but it should be
#   fairly comfortable & easy for spectroscopists to follow
function enact_init(O::op,ψ::psi)::Diagonal{Float64,Vector{Float64}}
   if O.a≠0
      out = O.v .* eh.(ψ.N).^O.a
   else
      out = fill(O.v,ψ.lng)
   end
   if O.b≠0; out .*= ns_el.(ψ.J,ψ.S,O.b,ψ.N)::Vector{Float64} ; end
   if O.c≠0; out .*= eh(ψ.S)^O.c ; end
   if O.d≠0; out .*= ψ.K .^ O.d ; end
   return Diagonal(out)
end

function torop(a,b,nf,ms)::SparseMatrixCSC{Float64,Int}
   out = 0.5.* (ms[1+b:end] .- nf*b).^a
   out = spdiagm(b=>out,-b=>reverse(out))
   return dropzeros!(out)
end

function enact(O::op,ψ::psi)::SparseMatrixCSC{Float64, Int}
#   out = O.v*O.f[1](ψ, O.p[1])::Diagonal{Float64,Vector{Float64}} #<- dispatch
   out = enact_init(O,ψ)
   @inbounds for i in 1:length(O.p)
#      mul!(out, out, O.f[i](ψ,O.p[i])::SparseMatrixCSC{Float64, Int}) #<- dispatch
      out *= O.f[i](ψ,O.p[i])::SparseMatrixCSC{Float64, Int}
   end
   if !isdiag(out) 
      tplus!(out) 
   end
   return out #<- dispatch 
end
#This allows the basis set to be distributed among a list of added operators
function enact(O::Vector{op},ψ::psi)::SparseMatrixCSC{Float64, Int}
   out = O[1]*ψ
   @inbounds for i in 1:length(O)
      out += enact(O[i],ψ)
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

units = Dict("MHz"=>1.,"cm-1"=>29979.2458,"kHz"=>"1e-3","Hz"=>"1e-6",
   "mHz"=>"1e-9","GHz"=>"1e3","THz"=>"1e6")


function ntop_enforce(O::op,lnfs)
   if size(O.tp,2) < lnfs
      O.tp = hcat(O.tp,zeros(2,lnfs - size(O.tp,2)))
   end
   return O
end

#Get it? \squarert? It's a more stable variant of sqrt that defaults to zero when
#   given a negative input value
eh(x::Real)::Float64 = x*(x+1)
□rt(x::Real)::Float64 =√(x*(x>zero(x)))
fh(x::Real,y::Real)::Float64 = □rt((x-y)*(x+y+1))
ns_el(j,s,p,n)::Float64 = (0.5*eh2(j) - eh2(n) - eh2(s))^p

#The four nice angular momentum operators
function nz(ψ::psi,p::Int)::Diagonal{Float64}
   Diagonal(ψ.K .^p)
end
function nt(ψ::psi,p::Int)::Diagonal{Float64}
   spdiagm(fill(eh(ψ.N)^p, ψ.lng))
end
function np(ψ::psi,p::Int)::SparseMatrixCSC{Float64, Int}
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



