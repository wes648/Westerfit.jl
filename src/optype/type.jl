
#lines I know won't work have been tagged with #<-----------

# string to function conversion:
# fn = getfield(Main, Symbol(str))

struct OpFunc{T <: Number, S<:AbPsi}
   f::FunctionWrapper{SparseMatrixCSC{T,Int}, Tuple{S,Int,Int}} # function
   l::Int # power / rank
   q::Int # top / component 
   OpFunc(T::Type,S::Type,f::Function,l::Int,q=0) = new{T,S}(f,l,q)
   function OpFunc(f::Function,l::Int,q=0)
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
      out = sand(out,tvs[O.q].vecs[:,:, σ2ind(ψ.σs[O.q], O.q, ψ.nfs[O.q])]) # <---------- need to determine the sigma index
      out = torsetter!(ψ,O.q,out)
   else
      @warn "unexpected condition in evalulation of tor op"
   end
   return out
end
#eval_top(op::OpFunc,ψ::TTPsi)::SparseMatrixCSC{T,Int} where {T<:Number} = op.f(ψ,op.l,op.q)
#eval_rop(op::OpFunc,ψ::RPsi)::SparseMatrixCSC{Float64,Int} = op.f(ψ,op.l,op.q)
#eval_top(op::OpFunc,ψ::TPsi)::SparseMatrixCSC{Float64,Int} = op.f(ψ,op.l,op.q)
struct Op
   nam::String #this is a name field, mostly for debugging
   rf::Vector{OpFunc} # vec of rot ops
   tf::Vector{OpFunc} # vec of tor ops
   stg::Int # stage
   Op(nam="E",rf=Vector{OpFunc}[],tf=Vector{OpFunc}[],stg=0) = new(nam,rf,tf,stg)
   Op(O::Op) = new(O.nam,O.rf,O.tf,O.stg)
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

#"""
#This function initiates the Eigs structure by calling the ctrl object.
#   It is mostly based on the number of stages and the nfold vector.
#   If stages == 1, one stage mode is initiated
#      |JSNKm₁,...,mₙ⟩ in one giant matrix (not advised due to run time)
#   If stages == 2, two stage mode is initiated
#      all tops in first stage, |m₁,...,mₙ⟩
#      |JSNK⟩⊗|vₜσ⟩ in the second stage
#   If stages == 3, three stage mode is initiated
#      each top considered individually: |m₁⟩,...,|mₙ⟩
#      then a toptop stage: |vₜ₁σ₁,...vₜₙσₙ⟩
#      finally tsr: |JSNK⟩⊗|vₜσ⟩
#"""
#function eigs_init(ctrl::Controls,nfold::Vector{Int})::Eigs
#   @warn "you need to double check this function wes"
#   if ctrl.stages==1
#      l3 = dgen(ctrl.mcalc)^length(nfold)*dgen(ctrl.s)*dgen(ctrl.Jmax)
#      s3 = length(σgen(nfold))
#      eigholder = Eigs(nothing, nothing, zeros(l3,l3,s3))
#   elseif stages==2
#      l2 = dgen(ctrl.mcalc)^length(nfold)
#      s = length(σgen(nfold))
#      l3 = (ctrl.vtcalc+1)*dgen(ctrl.s)*dgen(ctrl.Jmax)
#      eigholder = Eigs(nothing, zeros(l2,l2,s), zeros(l3,l3,s))
#   elseif stages==3
#      l1 = dgen(ctrl.mcalc)^length(nfold)
#      s1 = σcount(nfold[1])
#      for i ∈ 2:length(nfold)
#         s1 += nfold[i]
#      end
#      l2 = degen(ctrl.vtcalc)
#      s = length(σgen(nfold))
#      l3 = (ctrl.vtcalc+1)*dgen(ctrl.s)*dgen(ctrl.Jmax)
#      eigholder = Eigs(zeros(l1,l1,s1), zeros(l2,l2,s), zeros(l3,l3,s))
#   else
#      @warn "stages = $stages is not defined. Going to crash soon.
#               Julia will say that it's because eigholder isn't defined"
#   end
#   return eigholder
#end
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
#function stage_size2(ctrl::Controls)::Int
#   rscount = dgen(j)*dgen(ctrl.S)
#   σs = σcount(ctrl.NFOLD)
#   if isone(ctrl.stages)
#      l = rscount*(dgen(ctrl.mcalc)^length(ctrl.NFOLD))
#   elseif isone(length(ctrl.NFOLD)) && ctrl.stages > 1
#      lm = rscount*dgen(ctrl.mmcalc)
#      li = rscount*(ctrl.vtcalc + 1)
#      lv = rscount*(ctrl.vtmax + 1)
#      return new([ SubEigs(li,lm,σs) ],
#                 nothing,
#                 SubEigs(lv, li, σs) )
#   elseif ctrl.stages==2
#      l1 = rscount* (dgen(ctrl.mmcalc)^length(ctrl.NFOLD))
#      #l2 = rscount*(ctrl.vtrunc + 1)
#      l2 = rscount*(ctrl.vtcalc + 1)
#      l3 = rscount*(ctrl.vtmax + 1)
#      return new(nothing,
#                 SubEigs(l2,l1,σs),
#                 SubEigs(l3,l2,σs) )
#   elseif ctrl.stages==3
#      l1 = rscount* dgen(ctrl.mmcalc)
#      l2 = rscount*(ctrl.vtrunc + 1)
#      l3 = rscount*(ctrl.vtcalc + 1)
#      l4 = rscount*(ctrl.vtmax + 1)
#      return new([SubEigs(l2,l1,σs) for i ∈ 1:length(ctrl.NFOLD)],
#                 SubEigs(l3,l2,σs) ,
#                 SubEigs(l4,l3,σs) )
#   else
#      @warn "This number of stages ($(ctrl.stages)) has not been defined!"
#      new(nothing,nothing,nothing)
#   end # if
#   return l
#end # function

"""
Determines if previous stages are needed based on if the wavefunction matrix
   both exists (set by the stage keyword) and is nonzero (defined from said previous stage)
"""
stage_allow(x)::Bool = !isnothing(x) && !iszero(x.vecs)

