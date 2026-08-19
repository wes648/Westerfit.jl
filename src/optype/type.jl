
@kwdef mutable struct Controls
   apology::Bool = true
   RUNmode::String = "ESF"
   stages::Int = 0
   Irrep::String = "Ir"
   assign::String = "ram36"
   NFOLD::Vector{Int} = [0] # vector of symmetry folds of rotors
   S::Float64 = 0.0 # float for spin value could maybe turn into int for 2s but eh
   Jmax::Float64 = 0.0 # maximum J value
   mcalc::Int = 10 # maximum |m| value for free rotor basis, basis size will be 2mmax+1
   vtcalc::Int = 8 # maximum vt state output by second diag stage & to be used in third. basis size will be vtmax+1
   vtmax::Int = 0 # maximum vt state output by final diagonalization stage. basis size will be vtmax+1
   νrange::Vector{Float64} = [0.2; 40.0]
   TK::Float64 = 8. # temperature in Kelvin to be used in simulation
   INTthres::Float64 = 1e-6
   ExactHess::Bool = true # exact hess or not
   λlm0::Float64 = 0.001
   turducken::Int = 1
   maxiter::Int = 60
   BOLD::Int = 0
   REJECT::Float64 = 10.0
   goal::Float64 = 1.0
   overwrite::Bool = true 
end

struct OpFunc{S<:AbPsi}
   f::FunctionWrapper{SparseMatrixCSC{NUMTYPE,Int}, Tuple{S,Int,Int}} # function
   l::Int # power / rank
   q::Int # top / component 
   OpFunc(S::Type,f::Function,l::Int,q=0) = new{S}(f,l,q)
   function OpFunc(f::Function,l::Int,q=0)
      new{Tuple(first(methods(f)).sig.parameters[2:end])[1]}(f,l,q)
   end
end
struct MuFunc{T <: Number, S<:AbPsi}
   f::FunctionWrapper{SparseMatrixCSC{T,Int}, Tuple{S,S,Int,Int}} # function
   l::Int # power / rank
   q::Int # top / component 
   MuFunc(T::Type,S::Type,f::Function,l::Int,q=0) = new{T,S}(f,l,q)
   function MuFunc(f::Function,l::Int,q=0)
      new{Float64, Tuple(first(methods(f)).sig.parameters[2:end])[1]}(f,l,q)
   end
end

"""
eval_rop evaluates a rotational operator!
The arguments are op::OpFunc and ψ::RPsi
so function that acts on the rotational wavefunctions and said wavefunction
"""
eval_rop(op::OpFunc,ψ::RPsi)::SparseMatrixCSC{T,Int} where {T<:Number} = op.f(ψ,op.l,op.q)
#eval_top(op::OpFunc,ψ::TPsi)::SparseMatrixCSC{T,Int} where {T<:Number} = op.f(ψ,op.l,op.q)
"""
eval_tor evaluates a torsional operator!
The arguments are O::OpFunc, ψ::TTPsi, and tvs::Union{nothing,Matrix}
O acts on the torsional wavefunctions, ψ is the top-top wave function object
and tvs are the collection of one top wavefunctions or it's nothing.

No indexing can occur on tvs before the function is called as it might be nothing
"""
eval_top(O::OpFunc, ψ::TPsi)::SparseMatrixCSC{T,Int} where {T<:Number} = O.f(ψ, O.l, O.q)
function eval_top(O::OpFunc, ψ::TTPsi, tvs)#::SparseMatrixCSC{T,Int} where {T<:Number}
   out = O.f(ψ.tps[O.q], O.l, O.q)
   if isone(length(ψ.nfs)) # cases 0,1,2
      #nothing
   elseif length(ψ.nfs)>1 && isnothing(tvs) # cases 3,4,5
      out = torsetter!(ψ,O.q,out)
   elseif length(ψ.nfs)>1 && !isnothing(tvs) && iszero(tvs[O.q].vecs) # case 6
      #nothing
   elseif length(ψ.nfs)>1 && !isnothing(tvs) && !iszero(tvs[O.q].vecs) # case 7,8
      out = sand(out,tvs[O.q].vecs[:,:, σ2ind(ψ.σs[O.q], O.q, ψ.nfs[O.q])]) #
      out = torsetter!(ψ,O.q,out)
   else
      @warn "unexpected condition in evalulation of tor op"
   end
   return out
end
eval_rmu(ψb::RPsi,op::MuFunc,ψk::RPsi)::SparseMatrixCSC{T,Int} where {T<:Number} = op.f(ψb,ψk,op.l,op.q)
function eval_tmu(ψb::TTPsi, O::MuFunc, ψk::TTPsi, tvs)#::SparseMatrixCSC{T,Int} where {T<:Number}
   out = O.f(ψb.tps[O.q], ψk.tps[O.q], O.l, O.q)
   if isone(length(ψk.nfs)) # cases 0,1,2
      #nothing
   elseif length(ψk.nfs)>1 && isnothing(tvs) # cases 3,4,5
      out = torsetter!(ψk,O.q,out)
   elseif length(ψk.nfs)>1 && !isnothing(tvs) && iszero(tvs[O.q].vecs) # case 6
      #nothing
   elseif length(ψk.nfs)>1 && !isnothing(tvs) && !iszero(tvs[O.q].vecs) # case 7,8
      out = sand(out,tvs[O.q].vecs[:,:, σ2ind(ψk.σs[O.q], O.q, ψk.nfs[O.q])]) #
      out = torsetter!(ψk,O.q,out)
   else
      @warn "unexpected condition in evalulation of tor op"
   end
   return out
end

struct Op
   val::Float64
   rf::Vector{OpFunc{RPsi}}
   tf::Vector{OpFunc{TPsi}}
   Op(val=0.0,rf=Vector{OpFunc{RPsi}}[],tf=Vector{OpFunc{TPsi}}[]) = new(val,rf,tf)
   Op(val::Float64,rf::Vector{OpFunc{RPsi}},tf::Vector{OpFunc{TPsi}}) = new(val,rf,tf)
   Op(val::Float64,rf::Vector{OpFunc{RPsi}},tf::OpFunc{TPsi}) = new(val,rf,[tf])
   Op(val::Float64,rf::OpFunc{RPsi},tf::Vector{OpFunc{TPsi}}) = new(val,[rf],tf)
   Op(val::Float64,rf::OpFunc{RPsi},tf::OpFunc{TPsi}) = new(val,[rf],[tf])
   Op(O::Op) = new(O.val, O.rf, O.tf)
end
mutable struct Term
   const nam::String
   val::Float64
   const ops::Vector{Op}
   const scl::Float64
   const stg::Int
   const l::Int
   Term(nam::String,val::Float64,ops::Op,scl::Float64,stg::Int) = new(nam,val,[ops],scl,stg,1)
   Term(nam="E",val=0.0,ops=[Op()],scl=0.0,stg=0) = new(nam,val,ops,scl,stg,length(ops))
   Term(T::Term) = new(T.nam,T.ops,T.scl,T.stg,T.l)
end
struct Mu
   nam::String
   rf::Vector{MuFunc}
   tf::Vector{MuFunc}
   val::Float64
   Mu(nam="E",rf=Vector{MuFunc}[],tf=Vector{MuFunc}[],val=0.0) = new(nam,rf,tf,val)
   Mu(O::Mu) = new(O.nam,O.rf,O.tf)
end
#This is essentially a mutable version of the Eigen structure so I can rearrange & truncate as needed
# a is the number of saved eigen values
# b is the basis set size
# c is the number of σ states
mutable struct SubEigs
   vals::AbstractArray
   ders::Union{Nothing,AbstractArray}
   vecs::AbstractArray
   SubEigs(x::Eigen) = new(x.values,nothing,x.vectors)
   SubEigs(a::AbstractArray,c::AbstractArray) = new(a,nothing,c)
   SubEigs(a::AbstractArray,b::AbstractArray,c::AbstractArray) = new(a,b,c)
   """
   SubEigs can be initiated with 3 integers: number of states, basis size, and number of σs
      number of state should be less than or equal to the basis size
   """
   SubEigs(a::Int,b::Int,c::Int) = new(zeros(a,c), nothing, zeros(b,a,c))
end
mutable struct Eigs
   #the vector of subeigs refers to different tops
   # then vals[i,σ] is the ith energy of the σ symmetry
   # then vecs[:,i,σ] is the ith wavefunction of the σ symmetry
   top::Union{Nothing,Vector{SubEigs}}
   #vector inds will be sigma set
   # vals[i,σ] is the ith energy of the σ_set symmetry
   # vecs[:,i,σ] is the ith wavefunction of the σ_set symmetry
   ttp::Union{Nothing,SubEigs}
   #vib::Union{Nothing,SubEigs}
   #vector inds will be sigma set
   # vals[i,σ] is the ith energy of the σ_set symmetry
   # vecs[:,i,σ] is the ith wavefunction of the σ_set symmetry
   rst::Union{Nothing,SubEigs}
   function Eigs(ctrl::Controls)::Eigs
      rscount = sum(Int,2 .* (0.5*isodd(2ctrl.Jmax):ctrl.Jmax) .+ 1)*dgen(ctrl.S)
      rsd = dgen(ctrl.S)*dgen(ctrl.Jmax)
      σs = σcount(ctrl.NFOLD)
      if isone(ctrl.stages)
         lm = rsd*(dgen(ctrl.mcalc)^length(ctrl.NFOLD))
         lv = rscount*(ctrl.vtmax + 1) 
         return new( nothing, nothing, SubEigs(lv,lm,σs) )
#      elseif isone(length(ctrl.NFOLD)) && ctrl.stages > 1
#         lm = dgen(ctrl.mcalc)
#         li = (ctrl.vtcalc + 1)
#         lv = rscount*(ctrl.vtmax + 1)
#         return new([ SubEigs(li,lm,σs) ],
#                    nothing,
#                    SubEigs(lv, rsd*li, σs) )
      elseif ctrl.stages==2
         l1 = (dgen(ctrl.mcalc)^length(ctrl.NFOLD))
         l2 = (ctrl.vtcalc + 1)
         l3 = rscount*(ctrl.vtmax + 1)
         return new(nothing,
                    SubEigs(l2,l1,σs),
                    SubEigs(l3,rsd*l2,σs) )
      elseif ctrl.stages==3
         l1 =  dgen(ctrl.mcalc)
         l2 = (ctrl.vtcalc + 1)
         l3 = (ctrl.vtcalc + 1)
         l4 = rscount*(ctrl.vtmax + 1)
         return new([SubEigs(l2,l1,nth_σcount(ctrl.NFOLD[i],i)) for i ∈ 1:length(ctrl.NFOLD)],
                    SubEigs(l3,l2^length(ctrl.NFOLD),σs) ,
                    SubEigs(l4,rsd*l3,σs) )
      else
         @warn "This number of stages ($(ctrl.stages)) has not been defined!"
         new(nothing,nothing,nothing)
      end # if
   end # function
end # struct

"""
This function initiates the sparse zeros matrix for the new stage.
   Unfortunately it is hard coded.
"""
function stage_size(stage,stages,wvs)::Int
   println("stage is $stage")
   if stage==0# && stages ≥ 1
      return size(wvs.rst.vals,1)
   elseif (stage==2 && stages > 2) || (stage==1 && stages==2)
      return size(wvs.ttp.vals,1)
   elseif stage==1 && stages==3
      return size(wvs.top.vals,1)
   else
      @warn "stage = $stage is not defined. Going to crash soon."

   end
end
"""
Determines if previous stages are needed based on if the wavefunction matrix
   both exists (set by the stage keyword) and is nonzero (defined from said previous stage)
"""
stage_allow(x)::Bool = !isnothing(x) && !iszero(x.vecs)

