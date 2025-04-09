
import Base.*
import Base.^
import Base.+
import Base.-

Δlist(J::Real,S::Real)::Array{Int} = collect(Int(abs(J-S)):Int(J+S))
function ngen(J::Real,S::Real)::Vector{Int}
   ns = Δlist(J,S)
   out = fill(ns[1],2*ns[1]+1)
   for i in 2:length(ns)
      out = vcat(out,fill(ns[i],2*ns[i]+1))
   end
   return out
end
function kgen(J::Real,S::Real)::Vector{Int}
   ns = Δlist(J,S)
   out = collect(-ns[1]:ns[1])
   for i in 2:length(ns)
      n = ns[i]
      out = vcat(out,collect(-n:n))
   end
   return out
end
function msgen(nfold::Int,mcalc::Int,σ::Int)::StepRange{Int,Int}
   if iszero(nfold)
      marray = 0:0
   elseif isodd(nfold)
      lim = mcalc*nfold
      marray = (-lim+σ):nfold:(lim+σ)
   else #iseven(nfold)
      lim = floor(Int,lim/2)
      marray = -lim:floor(Int,nfold/2):lim
      marray = (-lim+σ):nfold:(lim+σ)
   end
   return marray
end
function msgen_indef(nf::Array{Int},mcalc::Int,σs::Array{Int})::Vector{StepRange{Int,Int}}
   out = Vector{StepRange{Int,Int}}(undef,length(nf))
   for j in eachindex(nf)
      out[j] = msgen(nf[j],mcalc,σs[j])
   end
   return out
end

Δlist2(J::Real,S::Real)::UnitRange{Int} = Int(abs(J-S)):Int(J+S)
kgen(ns::UnitRange{Int})::Vector{UnitRange{Int}} = [-n:n for n ∈ ns]

#tplus!(a::Diagonal)::SparseMatrixCSC{Float64, Int} = sparse(a)
#tplus!(a::Array{Float64,2})::Array{Float64,2} = hermitianpart!(2a)
function tplus!(a::SparseMatrixCSC{Float64, Int})::SparseMatrixCSC{Float64, Int}
   a .+= permutedims(a)
end

struct Psi
   J::Float64
   S::Float64
   N::UnitRange{Int}
   K::Vector{UnitRange{Int}}
   nf::Vector{Int}
   ms::StepRange{Int,Int}
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
         ms = msgen(nf,mc,σ)
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
         σ = [σ]
         nf = [nf]
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
         ms = msgen(nf,mc,σ)
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

struct OpFunc
   f::FunctionWrapper{SparseMatrixCSC{Float64,Int}, Tuple{Psi,Int}}
   p::Int
end
struct Op
   #this is a nam field
   nam::String
   #v for value as in parameter value
   #v::Float64
   #I've removed this field to make Ops fully immutable. the values will be in their own vector
   # since I need the mutability for the optimizer and have the hard coded 2nd order hamiltonian
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

################################################################################################
#   HI FUTURE WES!!! As of Apr7,'25 @ 13:43 You changed the definition of Op and need to fix 
#                      everything after this block
################################################################################################

function tarraysum(a::Array{Int,2},b::Array{Int,2})::Array{Int,2}
   if size(a,2) < size(b,2)
      a = hcat(a,zeros(Int,2,size(b,2)-size(a,2)))
   elseif size(b,2) < size(a,2)
      b = hcat(b,zeros(Int,2,size(a,2)-size(b,2)))
   end
   return a .+ b
end

function ur(n::Int)::SparseMatrixCSC{Float64, Int}
   out = Diagonal(vcat(fill(-√.5,n), 1.0, fill(√.5,n)))
   out += rotl90(Diagonal(vcat(fill(√.5,n), 0.0, fill(√.5,n))))
   return sparse(out)
end

function enact_init(O::Op,ψ::Psi,val)::SparseMatrixCSC{Float64,Int}#::Diagonal{Float64,Vector{Float64}}
   #out = O.a≠0 ? O.v .* eh.(ψ.N).^O.a : fill(O.v,ψ.lng)
   out = O.a≠0 ? n2(0.5*val, ψ.N, O.a) : fill(0.5*val,ψ.lng)
   if O.b≠0; out .*= ns_el2(ψ.J,ψ.S,ψ.N,O.b); end
   if O.c≠0; out .*= eh(ψ.S)^O.c ; end
   if O.d≠0; out .*= nz(ψ.K, O.d) ; end
   return spdiagm(out)
end

function torop(a::Int,b::Int,nf::Int,ms::StepRange{Int,Int})::SparseMatrixCSC{Float64,Int}
#performance of this function is deeply lacking
   if !iszero(b)
      @inbounds out = 0.5.* (ms[1+b:end] .- nf*b).^a
      out = spdiagm(b=>out)
      out[diagind(out,-b)] .= reverse!(-out)
      dropzeros!(out)
   elseif !iszero(a) && iszero(b)
      out = spdiagm(0=>ms .^a)
   else
      out = spdiagm(ones(size(ms)))
   end
   #@show out
   return out
end
function enact_tor(tp::Array{Int,2},ψ::Psi,
                   ut::SparseMatrixCSC{Float64,Int})::SparseMatrixCSC{Float64,Int}
   out = torop(tp[1,1],tp[2,1],ψ.nf[1],ψ.ms[1])
   if iszero(ψ.σ[1])
      out = dropzeros!(ut*out*ut)
   end
   @inbounds for i in 2:size(tp,2)
      part = torop(tp[1,i],tp[2,i],ψ.nf[i],ψ.ms[i])
      if iszero(ψ.σ[i])
         part = dropzeros!(ut*out*ut)
      end
      out = kron(part,out)
   end
   return out
end

function enact(O::Op,ψ::Psi,val::Float64,ur::SparseMatrixCSC{Float64,Int},
               ut::SparseMatrixCSC{Float64,Int})::SparseMatrixCSC{Float64,Int}
   out = enact_init(O,ψ,val)
   @inbounds for i in 1:length(O.rp)
      #part::SparseMatrixCSC{Float64,Int} = mateval(O.rf[i],ψ,O.rp[i])
      #out *= part
      #mul!(out,out,mateval(O.rf[i],ψ,O.rp[i]))
      #out *= O.rf[i](ψ::Psi,O.rp[i]::Int)::SparseMatrixCSC{Float64,Int} #<--------- CHECK
      out *= eval_op(O[i],ψ)
   end
   out = dropzeros!(ur*out*ur)
   #This block causes lots of JET errors
   if !iszero(O.tp) #O.tp ≠ zeros(Int,size(O.tp))
      part = enact_tor(O.tp,ψ)#<--------- CHECK
      out = kron(part,out)
   else
      out = kron(I(size(ψ.ms,1)),out)
   end
   #end said block of errors
   #if !isdiag(out) #<--------- CHECK
   #   tplus!(0.5*out) 
   #end
   return out #<- dispatch 
end
#This allows the basis set to be distributed among a list of added Operators
function enact(O::Vector{Op},ψ::Psi)::SparseMatrixCSC{Float64,Int}
   out = enact(O[1],ψ)
   @inbounds for i in 2:length(O)
      #part = enact(O[i],ψ)
      out += enact(O[i],ψ)::SparseMatrixCSC{Float64,Int}
   end
   return out
end
function torbuild(O::Vector{Op},ψ::Psi,stgs,msz,mc)::SparseMatrixCSC{Float64,Int}
   out = spzeros(msz,msz)
   @inbounds for i in 1:length(O)
      if isone(stgs[i])
         out += enact_tor(O[i].tp,ψ)*O[i].v
      elseif (stgs[i] < 0) && isone(stgs[i+stgs[i]])
         out += enact_tor(O[i].tp,ψ)*O[i+stgs[i]].v*O[i].v
      else
      end
   end
   if iszero(ψ.σ)
   U = ur(mc)
   mul!(out,U,out)
   mul!(out,out,U)
   end
   #@show out
   return dropzeros!(out)
end
function enact_stg2(O::Op,ψ::Psi,tvcs,mc)::SparseMatrixCSC{Float64,Int}
   #printstyled("start\n",color=:green)
   #@show O.v
   out = enact_init(O,ψ)
   @inbounds for i in 1:length(O.rp)
      out *= O.rf[i](ψ,O.rp[i])::SparseMatrixCSC{Float64,Int}
   end
   #@show out
   if !iszero(O.tp) #O.tp ≠ zeros(Int,size(O.tp))
      part = enact_tor(O.tp,ψ)
      if iszero(ψ.σ)
         U = ur(mc)
         mul!(part,U,part)
         mul!(part,part,U)
      end
      part = dropzeros!(sparse(tvcs' * part * tvcs))
      #@show part
      out = kron(part,out)
   else
      out = kron(I(size(tvcs,2)),out)
   end
   if !isdiag(out) 
      tplus!(0.5*out) 
   end
   #@show out
   #printstyled("stop\n",color=:red)
   return out #<- dispatch 
end
function h_stg2build!(Hmat,O::Vector{Op},ψ::Psi,stgs,siz,tvcs,mc)::SparseMatrixCSC{Float64,Int}
      @inbounds for i in eachindex(O) #1:length(O)
      if stgs[i] ≥ 2 #this is for future oddities 
         Hmat .+= enact_stg2(O[i],ψ,tvcs,mc)
      elseif stgs[i] < 0 && stgs[i] ≥ 2
         Hmat .+= enact_stg2(O[i],ψ,tvs,mc)*O[i+stgs[i]].v
      else
      end
   end
   return Hmat
end



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
=#
