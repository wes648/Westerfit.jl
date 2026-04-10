

#import Base.size
#size(x::Nothing,y::Int)::Int = 1

#tstage_check(ostg,stage,stages)::Bool = (ostg==stage) || ( 1  ≥ 1   &&   3  > 1)
ttstage_check(ostg,stage,stages)::Bool = (ostg==stage) || (ostg≥stage&&stages < 3)

function enact(O::Op,ψ::Psi,wvs::Eigs,val::Float64,UR::SparseMatrixCSC{Float64,Int},
               check=false)::SparseMatrixCSC#{T,Int} where T <: Number
#   out = enact_init(O,ψ.R,val) #potentially replace with just 0.5*val for improved readability
   #@show O.nam
   out = Diagonal(fill(0.5*val,ψ.R.lng))
   @inbounds for i ∈ eachindex(O.rf)
      out *= eval_rop(O.rf[i],ψ.R)
   end
   out = sand(out,UR)
   if !iszero(length(O.tf))#!isnothing(O.tf)
      if !isnothing(wvs.ttp) 
         tpart = sparse(I, size(wvs.ttp.vecs,1), size(wvs.ttp.vecs,1))
      else
         tpart = sparse(I,dgen(ψ.T.mc)^ψ.T.topcnt,dgen(ψ.T.mc)^ψ.T.topcnt)
      end
      @inbounds for i ∈ eachindex(O.tf)
         part = eval_top(O.tf[i], ψ.T, wvs.top)
         tpart = tpart*part
         if stage_allow(wvs.ttp)
#            println("hi!")
            tpart = sand(tpart, 
               wvs.ttp.vecs[:,:,nσfinder(O.tf[i].q, ψ.T.σs[O.tf[i].q], ψ.T.nfs)] )
         end #top-top if
      end #for i
      out = kron(tpart,out)
   elseif isnothing(wvs.ttp)
      out = kron(I(ψ.T.l), out)
   else # !isnothing(wvs.ttp.vals)
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
   @inbounds for i ∈ eachindex(O.tf)
      out *= eval_top(O.tf[i], ψ)
   end
   droptol!(out,2eps())
#   tplus!(out)
   return out
end 

function enact_tt(O::Op,ψ::TTPsi,wvs::Eigs,val::Float64)
   out = Diagonal(fill(0.5*val, size(wvs.ttp.vecs,1)))
   @inbounds for i ∈ eachindex(O.tf)
      tid = O.tf[i].q
      σid = σ2ind(ψ.tps[O.tf[i].q].σ, tid,ψ.tps[O.tf[i].q].nf) # <-------------------
      part = eval_top(O.tf[i], ψ, wvs.top) # <-----------------------
      out *= part
   end #for
   droptol!(out,2eps())
   #tplus!(out)
   return out
end

#processing 1 top hamiltonians
function one_topproc(wvs::Eigs,prms::Vector{Float64},ops,ϕ::TPsi,topid::Int,σid::Int,dmc::Int)::Eigs
   stage = 2#
   offset = hccount + 4*length(wvs.top)
# initialize Hamiltonian with size 2mc+1 x 2mc+1
# Hmat = spzeros(dmc,dmc)
   Hmat = htor2_hc(0.5 .*prms[hccount+4*(topid-1)+1:hccount+4*topid+1], ϕ)
   for i ∈ eachindex(ops)
      op = ops[i]
      if (op.stg == stage)&&(op.tf[1].q == topid)
         #ϕ = ψ.tps[op.tf[1].q]
         Hmat += enact_1t(op, ϕ, prms[i + offset])
      elseif (op.stg < 0)&&(op.tf[1].q == topid)&&(ops[i + op.stg].stg == stage)
         Hmat += enact_1t(op, ϕ, prms[i + offset]*prms[i + op.stg + offset])
      else
      end # stage if
   end # ops loop
   vals,vecs = diagwrap(tplus!(Hmat)) 
   l = size(wvs.top[topid].vals,1)
   wvs.top[topid].vals[:,σid] = vals[1:l]
   wvs.top[topid].vecs[:,:,σid] = vecs[:,1:l]
   return wvs
end
function stage_1tproc(wvs::Eigs,prms::Vector{Float64},ops,ctrl::Controls)::Eigs
   dmc = dgen(ctrl.mcalc)
   nfold = ctrl.NFOLD
   offset = hccount
   for i ∈ 1:length(nfold), j ∈ 1:nth_σcount(nfold[i],i) # <-------------
      σ = nth_σgen(nfold[i],i)[j]
      ϕ = TPsi(nfold[i],σ,ctrl.mcalc)
      wvs = one_topproc(wvs,prms,ops,ϕ,i,j,dmc)
   end # nested for
   return wvs
end

function stage_ttproc(wvs::Eigs,prms::Vector{Float64},ops::Vector{Op},
            ψ::TTPsi,σind::Int,ctrl::Controls)::Eigs
   #println("\nstart σ = $σind")
   stage = 1
   offset = hccount + 4*length(ctrl.NFOLD)
   l = ctrl.vtcalc+1
   #@show l
   if isnothing(wvs.top)
      ln = 2ctrl.vtcalc+1
      #Hmat = spzeros(size(wvs.ttp.vecs,1),size(wvs.ttp.vecs,1))
      Hmat = htor2_hc(0.5 .*prms[hccount+1:hccount+4], ψ.tps[1])
      for i ∈ 2:length(ctrl.NFOLD)
         Hmat = kron(sparse(I,ln,ln), Hmat ) + 
                kron( htor2_hc(prms[hccount+4i-3:hccount+4i], ψ.tps[i]),
                sparse(0.5I, ln^(i-1),ln^(i-1)),  ) # <-------------
      end
   else
      Hmat = 0.5 * wvs.top[1].vals[:,σ2ind(ψ,1)] 
      for i ∈ 2:length(ctrl.NFOLD)
         Hmat = kron(Hmat, ones(l) ) + 
                kron( fill(0.5, l^(i-1)), wvs.top[i].vals[:,σ2ind(ψ,i)]) # <-------------
      end
      #@show Hmat
      Hmat = spdiagm(Hmat)
   end
   for i ∈ eachindex(ops)
      if ttstage_check(ops[i].stg,stage,ctrl.stages)
         #part = enact_tt(ops[i],ψ,wvs,prms[i + offset])
         Hmat += enact_tt(ops[i],ψ,wvs,prms[i + offset])
      elseif ops[i].stg < 0 && ttstage_check(ops[i + ops[i].stg].stg,stage,ctrl.stages)
         Hmat += enact_tt(ops[i],ψ,wvs,prms[i + offset]*prms[i + op.stg + offset])
      end # stage if
   end # ops loop
   #if isone(σind)
   #   @show Hmat
   #end
   vals, vecs = diagwrap(tplus!(Hmat))
   wvs.ttp.vals[:,σind], wvs.ttp.vecs[:,:,σind] = vals[1:l], vecs[:,1:l]
   #println("end σ = $σind\n")
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
function stageproc0(ctrl,stage::Int,wvs::Eigs,prms::Vector{Float64},ops,ψ,σid)::Eigs
   offset = hccount + 4*length(ctrl.NFOLD)
   H = hrot2_hc(prms[1:4],ψ.R.N)
   if 0.0 < ψ.R.S < 1.0
      hsr!(H,prms[5:8], ψ.R)
   elseif 1.0 ≤ ψ.R.S
      hsr!(H,prms[5:8], ψ.R)
      hqu!(H,prms[9:11], ψ.R)
   else
      #nothing
   end
   U = ur(ψ.R.J,ψ.R.S)
   H = sand(H, U)
#   H += htsr2_hc(prms,wvs,ψ,σid)
   if !isnothing(wvs.ttp)
      H = kron(sparse(I(size(wvs.ttp.vals,1))), H )
      H[diagind(H)] .+= kron(wvs.ttp.vals[:,σid], fill(0.5,ψ.R.lng))
   elseif isone(ctrl.stages) && isone(length(ctrl.NFOLD)) && !iszero(ctrl.NFOLD[1])
      H = kron(sparse(I,ψ.T.l, ψ.T.l), H ) + kron(
         htor2_hc(prms[hccount+1:hccount+4], ψ.T.tps[1]), sparse(0.5I,ψ.R.lng,ψ.R.lng))
   elseif isone(ctrl.stages) && isone(length(ctrl.NFOLD)) && iszero(ctrl.NFOLD[1])
      # do nothing
   else
      @warn "what the fuck? you shouldn't be doing ntop 1 stage. this should have crashed"
   end
   for i ∈ eachindex(ops)
      #@show prms[i + offset]
      if ops[i].stg == stage || isone(ctrl.stages)
         H += enact(ops[i],ψ,wvs,prms[i + offset], U)
      elseif (ops[i].stg < 0)&& (ops[i + ops[i].stg].stg == stage || isone(ctrl.stages))
         H += enact(ops[i],ψ,wvs,prms[i + offset]*prms[i + ops[i].stg + offset], U)
      else
      end
   end
   tplus!(H)
   vals, vecs = diagwrap(H)
   if isone(ctrl.stages)&&isone(length(ctrl.NFOLD))
      perm = reshape(collect(1:ψ.R.lng).+ ψ.R.lng*(sortperm(by=abs,ψ.T.tps[1].ms) .-1)', ψ.R.lng*ψ.T.l)
      vecs .= vecs[perm,:]
   end
   perm = ram36_2stg_assign(vecs,ψ.R.J,ψ.R.S, ctrl.vtcalc, ctrl.vtmax)
   vldst = jinds(ψ.R.J, ψ.R.S, ctrl.vtmax+1) 
   vcdst = jinds(ψ.R.J, ψ.R.S, ctrl.vtcalc+1)
   wvs.rst.vals[vldst,σid], wvs.rst.vecs[1:size(vecs,1),vldst,σid] = vals[perm], vecs[:,perm]
   return wvs
end

function H_calc(ctrl::Controls,wvs::Eigs,prm_init,ops,jsσs::Matrix{Int})::Eigs
   prm = prm_proc(prm_init,length(ctrl.NFOLD))
   #@show prm
   σlist = unique(jsσs[:,2])
   σs = σgen(ctrl.NFOLD)
   if ctrl.stages==3
      #one top
      wvs = stage_1tproc(wvs,prm,ops,ctrl)
      println("one top stage done!")
   end
   if ctrl.stages ≥ 2
      # top-top
      for i ∈ σlist #1:size(σs,2)
         #println("now doing σ = $(σs[:,i])")
         if ctrl.stages > 2
            ψ = TTPsi(ctrl.NFOLD,σs[:,i],ctrl.mcalc,ctrl.vtcalc)
         else
            ψ = TTPsi(ctrl.NFOLD,σs[:,i],ctrl.mcalc)
         end
         wvs = stage_ttproc(wvs,prm,ops,ψ,i,ctrl)
         #@show wvs.ttp.vals
      end
      println("top-top stage done!")
   end
   for i ∈ 1:size(jsσs,1)
      σind = jsσs[i,2]
      j = 0.5*jsσs[i,1]
      ψ = Psi( RPsi(j, ctrl.S), TTPsi(ctrl.NFOLD,σs[:,σind],ctrl.mcalc) )
      wvs = stageproc0(ctrl,0,wvs,prm,ops,ψ,σind)
   end # σ loop
#   sparsify!(wvs.rst.vecs)
   return wvs
end

