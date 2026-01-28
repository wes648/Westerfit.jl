
import Base.*
import Base.^
import Base.+
import Base.-
import Base: convert

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
#This needs to get broken into 1 top and top-top structs
struct TPsi  <: AbPsi
   nf::Vector{Int}
   ms::Vector{StepRange{Int,Int}}
   σ::Vector{Int}
   lng::Int
   lb::Int
   function TPsi(nf::Vector{Int},σ::Vector{Int},mc=3)
      if length(nf) > 1
         ms = msgen(nf,mc,σ)
      else
         ms = [msgen(nf,mc,σ[1])]
      #   σ = [σ]
         nf = nf
      end
      lng = sum(length.(ms))
      #lb = 2mc+1
      new(nf, ms, σ, lng, 2mc+1)
   end
   function TPsi(nf::Vector{Int},σ::Int,mc=3)
      if isone(length(nf))
         ms = [msgen(nf[1],mc,σ)]
      end
      lng = sum(length.(ms))
      new(nf, ms, [σ], lng, 2mc+1)
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
   tops::Vector{TPsi_new}
   l::Int
   function TTPSi(top1::TPsi_new)
      new(1,[top1],top1.l)
   end
   function TTPsi(tops::TPsi_new...)
      new(length(tops),reduce(vcat,tops),reduce(+,tops.l))
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

struct OpFunc{T <: Number,S<:AbPsi}
   f::FunctionWrapper{SparseMatrixCSC{T,Int}, Tuple{S,Int}}
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
eval_rop(op::OpFunc,ψ::RPsi)::SparseMatrixCSC{T,Int} where {T<:Number} = op.f(ψ,op.p)
eval_top(op::OpFunc,ψ::TPsi)::SparseMatrixCSC{T,Int} where {T<:Number} = op.f(ψ,op.p)
struct Op
   #this is a name field, mostly for debugging
   nam::String
   #rf for rotational functions
   rf::Vector{OpFunc}
   #tf for torsional functions
   tf::Vector{OpFunc}
   #a,b,c,d are powers for N^2a (NS)^b S^2c Nz^d
   #this is to initialize the purely diagonal part of the matrix first
   #I can actually remove this entirely
   a::Int
   b::Int
   c::Int
   d::Int
   #forces the vector structure for even a single function
   function Op(nam::String,rf=Vector{OpFunc}[],tf=Vector{OpFunc}[],a=0,b=0,c=0,d=0)
      return new(nam,rf,tp,a,b,c,d)
   end
   function Op(nam::String,rf::Vector{OpFunc},abcd::Vector{Int},tf=Vector{OpFunc}[])
      return new(nam,rf,tp,abcd[1],abcd[2],abcd[3],abcd[4])
   end
   Op(O::Op) = new(O.nam,O.rf,O.tp,O.a,O.b,O.c,O.d)
end

mutable struct SubEigs
   vals::AbstractArray
   vecs::AbstractArray
end
mutable struct Eigs
   top::Union{Nothing,Vector{SubEigs}}
   ttp::Union{Nothing,SubEigs}
  #vib::Union{Nothing,SubEigs}
   rst::Union{Nothing,SubEigs}
end

#the below are for drafting and currently will not run
#lines I know won't work have been tagged with #<-----------
#I'm not super sure about the structure of a LT matrix ⊗ Dense
#I think the lower triangle would be fine?
function enact(O::Op,ψ::Psi,wvs::Eigs,val::Float64)::SparseMatrixCSC#{T,Int} where T <: Number
   out = enact_init(O,ψ.R,val) #potentially replace with just 0.5*val for improved readability
   @inbounds for i ∈ 1:length(O.rf)
      out *= eval_rop(O.rf[i],ψ.R)
   end
   if !isnothing(O.tf)
      tpart = sparse(1,1,1.0)
      @inbounds for i ∈ 1:length(O.tf)
         part = eval_top(O.tp,ψ.T)
         if !isnothing(wvs.top)
            part = sandwich(part, wvs.top[O[i].tf.q[σINDEX]]) #<-----------
         end
         torsetter!(out,XX,YY) #<-----------
         tpart = tpart*part
      end#for
      if !isnothing(wvs.ttp)
         tpart = sandwich(tpart, wvs.ttp[O[i].tf.q[XX]] )#<-----------
      end#top-top if
   else
      out = kron(I(ψ.T.lng),out)
   end#tor if
   #if
   #end vib
   dropzeros!(out)
   tplus!(out)#val gets multiplied by 0.5 in advance of this
   return out
end

function diagwrap(matrix)
   if eltype(H)<:Real #isreal(eltype(H))
      vals,vecs = eigen!(Symmetric(Matrix(H), :L))
   else
      vals,vecs = eigen!(Hermitian(Matrix(H), :L))
   end
   return vals, vecs
end

function stageproc(stage,wvs,prms,ops,ψ,stages)
   out = spzeros(????)#<-----------
   for i ∈ 1:length(ops)
      if stages[i]==stage
         out += eval(op,ψ,wvs,prms[i])
      elseif stages[i] < 0
         out += eval(op,ψ,wvs,prms[i])*prms[i-stages[i]]
      else
      end
   end
   if iszero(stage)
      U = kron(sparse(1.0I,ψ.T.lng,ψ.T.lng), ur(ψ.R.J,ψ.R.S)) #<-----------
      H = droptol!(sand(H,U),2*eps())
      vals, vecs = diagwrap(out)
      ASSIGN #<-----------
   else# all other stages use energetic based assignments
      vals, vecs = diagwrap(out)
   end
   wvs.vals, wvs.vecs = vals, vecs
   return wvs
end

function H_calc()