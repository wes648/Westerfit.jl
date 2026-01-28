
abstract type AbPsi end

struct RPsi <: AbPsi
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
struct TPsi_new  <: AbPsi
   nf::Int
   ms::StepRange{Int,Int}
   σ::Int
   l::Int
   function TPsi(nf::Int,σ::Int,mc=3)
      ms = msgen(nf,mc,σ)
      new(nf,ms,σ,length(ms))
   end
end
struct TTPsi 
   topcnt::Int
   tps::Vector{TPsi_new}
   nfs::Vector{Int}
   σs::Vector{Int}
   l::Int
   function TTPSi(top1::TPsi_new)
      new(1,[top1],top1.nf,top1.σ,top1.l)
   end
   function TTPsi(tops::TPsi_new...)
      new(length(tops),reduce(vcat,tops),reduce(vcat,tops.nf),reduce(vcat,tops.σ),reduce(+,tops.l))
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
   T::TTPsi
end

convert(T::Type{Psi},ϕ::RPsi) = Psi(ϕ,TPsi(0,0,0))
convert(T::Type{Psi},ϕ::TPsi) = Psi(RPsi(0),ϕ)
