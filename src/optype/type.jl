
#lines I know won't work have been tagged with #<-----------

@kwdef struct Controls
   NFOLD::Vector{Int} = [0]
   S::Float64 = 0.0
   TK::Float64 = 8.0
   mcalc::Int = 8
   mmax::Int = 6
   vtmax::Int = 0
   Jmax::Float64 = 0.
   apology::Bool = true
   νmin::Float64 = 0.0
   νmax::Float64 = 40.0
   INTthres::Float64 = 0.0001
   λlm0::Float64 = 0.001
   RUNmode::String = "ESF"
   turducken::Int = 1
   maxiter::Int = 60
   overwrite::Bool = true 
   assign::String = "ram36"
   REJECT::Float64 = 10.0
   Irrep::String = "Ir"
   goal::Float64 = 1.0
   stages::Int = 1 
   ctbk::Vector{Int} = zeros(Int,2)
   sobk::Vector{Int} = zeros(Int,2)
   inbk::Vector{Int} = zeros(Int,2)
   opbk::Vector{Int} = zeros(Int,2)
end

struct OpFunc{T <: Number, S<:AbPsi}
   f::FunctionWrapper{SparseMatrixCSC{T,Int}, Tuple{S,Int}} # function
   l::Int # power / rank
   q::Int # component / top
   OpFunc(T::Type,f::Function,l::Int,q=0) = new{T}(f,l,q)
   OpFunc(f::Function,l::Int,q=0) = new{Float64}(f,l,q)
end
eval_rop(op::OpFunc,ψ::RPsi)::SparseMatrixCSC{T,Int} where {T<:Number} = op.f(ψ,op.l,op.q)
eval_top(op::OpFunc,ψ::TPsi)::SparseMatrixCSC{T,Int} where {T<:Number} = op.f(ψ,op.l,op.q)
struct Op
   nam::String #this is a name field, mostly for debugging
   rf::Vector{OpFunc} # vec of rot ops
   tf::Vector{OpFunc} # vec of tor ops
   stg::Int # stage
   Op(nam::String,rf=Vector{OpFunc}[],tf=Vector{OpFunc}[],stg=0) = new(nam,rf,tf,stg)
   Op(O::Op) = new(O.nam,O.rf,O.tf,O.stg)
end

#This is essentially a mutable version of the Eigen structure so I can rearrange & truncate as needed
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

"""
This function initiates the Eigs structure by calling the ctrl object.
   It is mostly based on the number of stages and the nfold vector.
   If stages == 1, one stage mode is initiated
      |JSNKm₁,...,mₙ⟩ in one giant matrix (not advised due to run time)
   If stages == 2, two stage mode is initiated
      all tops in first stage, |m₁,...,mₙ⟩
      |JSNK⟩⊗|vₜσ⟩ in the second stage
   If stages == 3, three stage mode is initiated
      each top considered individually: |m₁⟩,...,|mₙ⟩
      then a toptop stage: |vₜ₁σ₁,...vₜₙσₙ⟩
      finally tsr: |JSNK⟩⊗|vₜσ⟩
"""
function eigs_init(ctrl::Controls,nfold::Vector{Int})::Eigs
   if ctrl.stages==1
      l3 = dgen(ctrl.mcalc)^length(nfold)*dgen(ctrl.s)*dgen(ctrl.jmax)
      s3 = length(σgen(nfold))
      eigholder = Eigs(nothing, nothing, zeros(l3,l3,s3))
   elseif stages==2
      l2 = dgen(ctrl.mcalc)^length(nfold)
      s = length(σgen(nfold))
      l3 = (ctrl.vtcalc+1)*dgen(ctrl.s)*dgen(ctrl.jmax)
      eigholder = Eigs(nothing, zeros(l2,l2,s), zeros(l3,l3,s))
   elseif stages==3
      l1 = dgen(ctrl.mcalc)^length(nfold)
      s1 = σcount(nfold[1])
      for i ∈ 2:length(nfold)
         s1 += nfold[i]
      end
      l2 = degen(vctrl.tcalc)
      s = length(σgen(nfold))
      l3 = (ctrl.vtcalc+1)*dgen(ctrl.s)*dgen(ctrl.jmax)
      eigholder = Eigs(zeros(l1,l1,s1), zeros(l2,l2,s), zeros(l3,l3,s))
   else
      @warn "stages = $stages is not defined. Going to crash soon.
               Julia will say that it's because eigholder isn't defined"
   end
   return eigholder
end
"""
This function initiates the sparse zeros matrix for the new stage.
   Unfortunately it is hard coded.
"""
function stage_size(stage,stages,wvs)::Int
   if stage==0# && stages ≥ 1
      return size(wvs.rst,1)
   elseif (stage==2 && stages > 2) || (stage==1 && stages==2)
      return size(wvs.ttp,1)
   elseif stage==1 && stages==3
      return size(wvs.top,1)
   else
      @warn "stage = $stage is not defined. Going to crash soon."

   end
end

"""
Determines if previous stages are needed based on if the wavefunction matrix
   both exists (set by the stage keyword) and is nonzero (defined from said previous stage)
"""
stage_allow(x)::Bool = !isnothing(x) && !iszero(x)

