
function enact(O::Op,lng::Int,ψ::Psi,wvs::Eigs,val::Float64)::SparseMatrixCSC#{T,Int} where T <: Number
#   out = enact_init(O,ψ.R,val) #potentially replace with just 0.5*val for improved readability
   out = Diagonal(fill(0.5*val,lng))
   @inbounds for i ∈ 1:length(O.rf)
      out *= eval_rop(O.rf[i],ψ.R)
   end
   out = sand(out,ur(ψ.R.J, ψ.R.S))
   if !isnothing(O.tf)
      tpart = sparse(1,1,1.0)
      @inbounds for i ∈ 1:length(O.tf)
         part = eval_top(O.tf[i],ψ.T)
         if stage_allow(wvs.top) 
            loc = nσfinder(O.tf[i].q, ψ.T.σs[top_id])
            part = sand(part, wvs.top[loc].vecs)
         end
         torsetter!(ψ.T, O.tf[i].q, part)
         tpart = tpart*part
      end#for
      if stage_allow(wvs.ttp)
         tpart = sand(tpart, wvs.ttp[σfinder(O[i].tf.q, ψ.T.σs[O.tf[i].q], ψ.T.nfs)] )
      end#top-top if
      out = kron(tpart,out)
   else
      out = kron(I(ψ.T.lng),out)
   end#tor if
   #if
   #end vib
   droptol!(out,2eps())
   tplus!(out)#val gets multiplied by 0.5 in advance of this
   return out
end

function enact_1t(Hmat,O::Op,nfold::Vector{Int},lng::Int,topid::Int,ψ::TTPsi,val::Float64)
   out = Diagonal(fill(0.5*val,lng))
   @inbounds for i ∈ 1:length(O.tf[i])
      out *= eval_top(O.tf[i], ψ)
   end
   droptol!(out,2eps())
   tplus!(out)
   return out
end

function enact_tt()
end


"""
stageproc is the stage processor. It build the ℋ matrix for a specific stage. 
inputs are: stage id #
            wvs is the wavefunctions (eigs structure)
            prms is the parameter vector
            ops is vector of operators
            ψ is the quantum number structure
            stages is the vector if stages for the operators
Outputs the Eig structure
"""
function stageproc(stage::Int,wvs::Eigs,prms::Vector{Float64},ops,ψ,stages)::Eigs
   l = stage_size(stage,stages,wvs)
   H = spzeros(l,l)
   for i ∈ 1:length(ops)
      if ops[i].stg == stage
         H += eval(ops[i],l,ψ,wvs,prms[i])
      elseif (ops[i].stg < 0)&& (ops[i + ops[i].stg].st == stage)
         H += eval(ops[i],ψ,wvs,prms[i])*prms[i + ops[i].stg ]
      else
      end
   end
   if iszero(stage)
      U = kron(sparse(1.0I,ψ.T.lng,ψ.T.lng), ur(ψ.R.J,ψ.R.S))
      H = droptol!(sand(H,U),2*eps())
      vals, vecs = diagwrap(H)
      ASSIGN #<-----------
   else# all other stages use energetic based assignments
      vals, vecs = diagwrap(H)
   end
   wvs.vals, wvs.vecs = vals, vecs
   return wvs
end

function stage_1tproc()::Eigs
   Hmat = spzeros(dgen(mcalc),dgen(mcalc),length(nfold))
   for i ∈ 1:length(ops)
      Hmat[:,:,nσfinder(topid, ψ.σ, nfold)] += enact_1t()
   end
   for i ∈ 1:length(nfold)
      vals, vecs = diagwrap(Hmat[:,:,i])
      wvs.top[i].vals = vals[1:ctrl.vtcalc]
      wvs.top[i].vecs = vecs[:,:, 1:ctrl.vtcalc]
   end
   return wvs
end


function H_calc(ctrl::Controls,wvs::Eigs,vals)::Eigs
   for  stage ∈ 1:ctrl.stages
      wvs = stageproc(stage,wvs,prms,ops,ψ,ctrl.stages)
   end
   wvs
end

#=
if stages==1
   eval all
elseif stages==2
   eval !iszero(stage) ops
   eval iszero(stage) ops
=#