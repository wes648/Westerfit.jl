
ttstage_check(ostg,stage,stages)::Bool = (ostg==stage) || (ostg≥stage&&stages>stage)

function enact(O::Op,ψ::Psi,wvs::Eigs,val::Float64,check=false)::SparseMatrixCSC#{T,Int} where T <: Number
#   out = enact_init(O,ψ.R,val) #potentially replace with just 0.5*val for improved readability
   out = Diagonal(fill(0.5*val,ψ.R.lng))
   @inbounds for i ∈ eachindex(O.rf)
      out *= eval_rop(O.rf[i],ψ.R)
   end
   out = sand(out,ur(ψ.R.J, ψ.R.S))
   if !iszero(length(O.tf))#!isnothing(O.tf)
      tpart = sparse(I,ψ.T.l,ψ.T.l)
      @inbounds for i ∈ eachindex(O.tf)
         part = eval_top(O.tf[i], ψ.T.tps[O.tf[i].q])
         if stage_allow(wvs.top) 
            loc = nσfinder(O.tf[i].q, ψ.T.σs[top_id])
            part = sand(part, wvs.top[loc].vecs)
         end
         torsetter!(ψ.T, O.tf[i].q, part)
         tpart = tpart*part
         if stage_allow(wvs.ttp)
            tpart = sand(tpart, 
               wvs.ttp.vecs[:,:,nσfinder(O.tf[i].q, ψ.T.σs[O.tf[i].q], ψ.T.nfs)] )
         end #top-top if
      end #for i
      out = kron(tpart,out)
   else
      out = kron(I(size(wvs.ttp.vals,1)),out)
   end #tor if
   #if
   #end vib
   droptol!(out,5eps())
   if check
      tplus!(out)
   end #check if
#   tplus!(out)#val gets multiplied by 0.5 in advance of this
   return out
end

function enact_1t(O::Op,ψ::TPsi,val::Float64)
   out = Diagonal(fill(0.5*val,ψ.l))
   @inbounds for i ∈ eachindex(O.tf[i])
      out *= eval_top(O.tf[i], ψ, nothing)
   end
   droptol!(out,2eps())
#   tplus!(out)
   return out
end 

function enact_tt(O::Op,ψ::TTPsi,wvs::Eigs,val::Float64)
   out = Diagonal(fill(0.5*val,ψ.l))
   @inbounds for i ∈ eachindex(O.tf)
      topid = O.tf[i].q
      part = eval_top(O.tf[i], ψ.tps[topid], )
      if stage_allow(wvs.top) 
         loc = nσfinder(O.tf[i].q, ψ.tps.σs[top_id])
         part = sand(part, wvs.top[loc].vecs)
      end
      torsetter!(ψ, O.tf[i].q, part)
      out *= part
   end #for
   droptol!(out,2eps())
   tplus!(out)
   return out
end

#processing 1 top hamiltonians
function one_topproc(wvs::Eigs,prms::Vector{Float64},ops,ϕ::TPsi,topid::Int,σid::Int,dmc::Int)::Eigs
   stage = 2#
# initialize Hamiltonian with size 2mc+1 x 2mc+1
   Hmat = spzeros(dmc,dmc)
   for i ∈ eachindex(ops)
      op = op[i]
      if (op.stg == stage)&&(op.tf[1].q == topid)
         #ϕ = ψ.tps[op.tf[1].q]
         Hmat += enact_1t(op, ϕ, prms[i])
      elseif (op.stg < 0)&&(op.tf[1].q == topid)&&(ops[i + op.stg].st == stage)
         Hmat += enact_1t(op, ϕ, prms[i]*prms[i + op.stg])
      else
      end # stage if
   end # ops loop
   wvs.top[topid].vals[:,σid], wvs.top[topid].vecs[:,:,σid] = diagwrap(tplus!(Hmat))
   return wvs
end
function stage_1tproc(wvs::Eigs,ops,ctrl::Controls)::Eigs
   dmc = dgen(ctrl.mcalc)
   nfold = ctrl.NFOLD
   for i ∈ 1:length(nfold), j ∈ 1:nth_σcount(nfold[i],i) # <-------------
      σ = nth_σgen(nfold[i],i)
      ϕ = TPsi(nfold[i],σ,ctrl.mcalc)
      wvs = one_topproc(wvs,prms,ops,ϕ,i,j,dmc)
   end # nested for
   return wvs
end

function stage_ttproc(wvs::Eigs,prms::Vector{Float64},ops,ψ,σind::Int,ctrl)::Eigs
   stage = 1
   Hmat = spzeros(ψ.l,ψ.l)
   for i ∈ eachindex(ops)
      if ttstage_check(ops[i].stg,stage,ctrl.stages)
         Hmat += enact_tt(ops[i],ψ,wvs,prms[i])
      elseif ops[i].stg < 0 && ttstage_check(ops[i + ops[i].stg].stg,stage,ctrl.stages)
         Hmat += enact_tt(ops[i],ψ,wvs,prms[i]*prms[i + op.stg])
      end # stage if
   end # ops loop
   #@show Hmat
   vals, vecs = diagwrap(tplus!(Hmat))
   #@show vals
   wvs.ttp.vals[:,σind], wvs.ttp.vecs[:,:,σind] = vals[1:ctrl.vtcalc+1], vecs[:,1:ctrl.vtcalc+1]
   return wvs
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
function stageproc0(ctrl,stage::Int,wvs::Eigs,prms::Vector{Float64},ops,ψ,stages,σid)::Eigs
   if !isnothing(wvs.ttp)
      H = spdiagm(kron(wvs.ttp.vals[:,σid],ones(ψ.R.lng)))
   else
      #l = ψ.R.lng*ψ.T.l
      H = spzeros(ψ.R.lng*ψ.T.l, ψ.R.lng*ψ.T.l)
   end
   for i ∈ eachindex(ops)
      if ops[i].stg == stage
         H += enact(ops[i],ψ,wvs,prms[i])
#         part = enact(ops[i],ψ,wvs,prms[i]); @show size(part)
#         H += part 
      elseif (ops[i].stg < 0)&& (ops[i + ops[i].stg].st == stage)
         H += enact(ops[i],ψ,wvs,prms[i]*prms[i + ops[i].stg ])
      else
      end
   end
   tplus!(H)
   #if assigncheck # energy caculation!
      vals, vecs = diagwrap(H)
      if isone(ctrl.stages)&&isone(length(ctrl.NFOLD))
         perm = reshape(collect(1:ψ.R.lng).+ ψ.R.lng*(sortperm(by=abs,ψ.T.tps[1].ms) .-1)', ψ.R.lng*ψ.T.l)
         vecs .= vecs[perm,:]
      end
      perm = ram36_2stg_assign(vecs,ψ.R.J,ψ.R.S, ctrl.vtcalc, ctrl.vtmax)
      vldst = jinds(ψ.R.J,ψ.R.S,ctrl.vtmax+1) 
      vcdst = jinds(ψ.R.J,ψ.R.S,ctrl.vtcalc+1)
      wvs.rst.vals[vldst,σid], wvs.rst.vecs[1:size(vecs,1),vldst,σid] = vals[perm], vecs[:,perm]
   #else # derivatives!
   #   wvs.rst.der[:,σid] = 
   #end
   return wvs
end

function H_calc(ctrl::Controls,wvs::Eigs,prm,ops,stages)::Eigs
   σs = σgen(ctrl.NFOLD)
   if ctrl.stages==3
      #one top
      wvs = stage_1tproc(wvs,prm,ops,ctrl)
   end
   if ctrl.stages ≥ 2
      # top-top
      for i ∈ 1:size(σs,2)
         ψ = TTPsi(ctrl.NFOLD,σs[:,i],ctrl.mcalc)
         wvs = stage_ttproc(wvs,prm,ops,ψ,i,ctrl)
      end
   end
   for i ∈ 1:size(σs,2)
      ϕ = TTPsi(ctrl.NFOLD,σs[:,i],ctrl.mcalc)
      for j ∈ 0.5*isodd(2*ctrl.S):ctrl.Jmax
         ψ = Psi( RPsi(j, ctrl.S), ϕ )
         wvs = stageproc0(ctrl,0,wvs,prm,ops,ψ,stages,i)
      end # j loop
   end # σ loop
#   sparsify!(wvs.rst.vecs)
   return wvs
end

#=
if stages==1
   eval all
elseif stages==2
   eval !iszero(stage) ops
   eval iszero(stage) ops
=#
