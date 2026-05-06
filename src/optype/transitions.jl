
function enact_mu(O::MuOp,val::Float64,ψb::Psi,σb::Int,ψk::Psi,σk::Int,wvs::Eigs)
   out = val
   for i ∈ eachindex(O.rf)
      out *= eval_rmu(ψb.R,O.rf[i],ψk.R)
   end
   if !iszero(length(O.tf)) # how to eval the tor_op
      if !isnothing(wvs.ttp) # 2 & 3 stage mode
         tpart = sparse(I, size(wvs.ttp.vecs,1), size(wvs.ttp.vecs,1))
      else # one stage mode
         tpart = sparse(I,dgen(ψ.T.mc)^ψ.T.topcnt,dgen(ψ.T.mc)^ψ.T.topcnt)
      end # end wvs.ttp if
      for i ∈ eachindex(O.tf)
         part = eval_tmu(ψb.R,O.rf[i],ψk.R, wvs.top)
         tpart = tpart*part
         if stage_allow(wvs.ttp)
            # ⟨ ψ' | O_t | ψ ⟩
            tpart = 
               wvs.ttp.vecs[:,:,nσfinder(O.tf[i].q, ψb.T.σs[O.tf[i].q], ψb.T.nfs)]' *
               tpart *
               wvs.ttp.vecs[:,:,nσfinder(O.tf[i].q, ψk.T.σs[O.tf[i].q], ψk.T.nfs)] 
         end #top-top if
      end #for i
      out = kron(tpart,out)
   elseif isnothing(wvs.ttp) # one stage no op
      out = kron(I(ψ.T.l), out)
   else # !isnothing(wvs.ttp.vals) # two stage no op
      out = kron(wvs.ttp.vecs[:,:,σbid]' * wvs.ttp.vecs[:,:,σkid],out) # ⟨ ψt' | ψt ⟩
   end #tor if
   # ⟨ ψtr' | O_f | ψtr ⟩
   out = sparsify!(
      wvs.rst.vecs[1:size(out,2), jinds(ψb.R), σb]' *
      out * 
      wvs.rst.vecs[1:size(out,1), jinds(ψK.R), σk]
   ,1e-5)
   return out
end

function qtot_calc(jmax::Float64,s::Float64,TK::Float64,wvs::Eigs)::Float64
   jdlist = mapreduce(x->fill(2x+1,Int(2x)+1), append!, jlister(jmax,s))
   out = exp.(-wvs.rst.vals .* (kb*TK)^-1 )) .* jdlist
   if size(out,2) > 1
      out[:,2:end] .* 2
   end
   return (2s+1)*sum!(out)
end

function μ2_build(Ops::Vector{MuOp},vals::Vector{Float64},ψb::Psi,σb::Int,ψk::Psi,σk::Int,wvs::Eigs)
   μmat = spzeros( ????? ) # <----------
   for i ∈ eachindex(Ops)
      μmat += enact_mu(Ops[i], vals[i], ψb,σb,ψk,σk)
   end
   μmat .^= 2
   if σb==σk && ψb.R.J==ψk.R.J
      μmat = sparse( UpperTriangular(μmat) - Diagonal(μmat) )
   end
   return μmat
end

function μ_proc!(ints,Ops,vals,ψb,σb,ψjk,σk,wvs)
   binds = ???? # <----------
   kinds = ???? # <----------
   rind, cind, tints = findnz( μ2_build(Ops,vals, ψb,σb, ψk,σk, wvs) )
   for l ∈ eachindex(tints)
      ints[binds[rind[l]], kinds[cind[l]] ] = tints[l]
   end
end

function thermeff(ν,El,kT,Q)
   x = exp(-El/(kT*csl))
   return ν*x*(1.0-x) / Q*csl
end

function μσbσk_stage(ctrl, kT,Q, ϕb,σb, ϕk,σk, wvs)
   ints = spzeros(size(wvs.rst.vals,1), size(wvs.rst.vals,1))
   if σb==σk
      jbjk = ????? # <----------
   elseif σb≠σk
      jbjk = ????? # <----------
   else
      @warn "wtf intensity calc ⟨σ'|σ⟩ pairing issue"
   end
   for i ∈ 1:size(jbjk)
      ψb = Psi(RPsi(jbjk[i,1],ctrl.S), ϕb)
      ψk = Psi(RPsi(jbjk[i,2],ctrl.S), ϕk)
      μ_proc!(ints,Ops,vals,ψb,σb,ψk,σk,wvs)
   end
   outfs = zeros(Float64, nnz(ints),3)
   outis = zeros(Int, nnz(ints),4)
   rinds, cinds, = findnnz(ints)
   for i ∈ 1:nnz(ints)
      r = rinds[i]
      c = cinds[i]
      ν = wvs.rst.vals[r,σb] - wvs.rst.vals[c,σk]
      if ctrl.νrange[2] > ν > ctrl.νrange[1]
         outfs[i,1] = ν
         outfs[i,2] = ints[i]*thermeff(ν,wvs.rst.vals[c,σk], kt,Q)
         outfs[i,3] = wvs.rst.vals[c,σk]/csl
         outis[i,:] = [r σb c σk]
      elseif ctrl.νrange[2] > -ν > ctrl.νrange[1]
         outfs[i,1] = -ν
         outfs[i,2] = ints[i]*thermeff(-ν,wvs.rst.vals[r,σb], kt,Q)
         outfs[i,3] = wvs.rst.vals[r,σb]/csl
         outis[i,:] = [c σk r σb]
      else
      end
   end
   filt = findall(!iszero, outfs[:,1])
   outis = outis[filt,:] 
   outfs = outfs[filt,:]
   filt = findall(x -> x>ctrl.INTthres, outfs[:,2])
   return outis[filt,:], outfs[filt,:]
end

function tracalc(ctrl,wvs,molnam)
   inds = zeros(Int, 0,4)
   frqs = zeros(Float64, 0,3)
   Q = ????
   kT = ctrl.TK * kb
   for σb ∈ 1:σcount, σk ∈ 1:σb
      ϕb = ????? # <------
      ϕk = ????? # <------
      tinds, tfrqs = μσbσk_stage(ctrl, kT,Q, ϕb,σb, ϕk,σk, wvs)
      inds = vcat(inds,tinds)
      frqs = vcat(frqs,tfrqs)
   end
   perm = sortperm(frqs)
   frqs = frqs[perm,:]
   inds = inds[perm,:]
   writefreqs(frqs,inds)
   return frqs, inds
end