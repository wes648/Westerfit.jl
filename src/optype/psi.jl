
abstract type AbPsi end

struct RPsi <: AbPsi
   J::Float64
   S::Float64
   N::UnitRange{Int}
   K::Vector{UnitRange{Int}}
   lng::Int
   function RPsi(J::Number,S::Number)
      J = Float64(J)
      S = Float64(S)
      N = Δlist2(J,S)
      K = kgen(N)
      lng = Int((2J+1)*(2S+1))
      new(J,S,N,K,lng)
   end
   function RPsi(N::Int)
      J = Float64(N)
      N = Δlist2(J,0.0)
      K = kgen(N)
      lng = Int(2J+1)
      new(J,0.0,N,K,lng)
   end
end
struct TPsi <: AbPsi
   nf::Int
   ms::StepRange{Int,Int}
   σ::Int
   l::Int
   function TPsi(nf::Int,σ::Int,mc=3)
      ms = msgen(nf,mc,σ)
      new(nf,ms,σ,length(ms))
   end
end
struct TTPsi <: AbPsi
   topcnt::Int
   tps::Vector{TPsi}
   nfs::Vector{Int}
   σs::Vector{Int}
   mc::Int
   l::Int
   function TTPSi(top1::TPsi)
      new(1, [top1], top1.nf, top1.σ, floor(Int,0.5*(top1.l-1)), top1.l )
   end
   function TTPsi(tops::TPsi...)
      new(length(tops),reduce(vcat,tops),reduce(vcat,tops.nf),
         reduce(vcat,tops.σ),reduce(x -> floor(Int, 0.5*(x-1)), tops.l),
         reduce(+,tops.l))
   end
   function TTPsi(nfs::Vector{Int},σs::Vector{Int},mc::Int)
      topcnt = length(nfs)
      tps = Vector{TPsi}(undef,length(nfs))
      for i ∈ 1:topcnt
         tps[i] = TPsi(nfs[i],σs[i],mc)
      end
      new(topcnt, tps,nfs, σs, mc, dgen(mc)^topcnt)
   end
   TTPsi(nfs::Int,σs::Int,mc::Int) = TTPsi([nfs],[σs],mc)
end
struct VPsi
   #each column refer to a vibrational polyad state
   #each row refers to the quanta of a given mode
   #thus nus[a,b] is the quanta of mode a in state b
   nus::Matrix{Int}
end
mutable struct Psi
   R::Union{RPsi,Nothing}
   T::Union{TTPsi,Nothing}
end

convert(T::Type{Psi},ϕ::RPsi) = Psi(ϕ,TPsi(0,0,0))
convert(T::Type{Psi},ϕ::TPsi) = Psi(RPsi(0),ϕ)