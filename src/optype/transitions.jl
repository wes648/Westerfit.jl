function qtot_calc(jmax::Float64,s::Float64,TK::Float64,wvs::Eigs)::Float64
   jdlist = mapreduce(x->fill(x+1,Int(x)+1), append!, jlister(jmax,s))
   out = exp.(-wvs.rst.vals .* (kb*TK)^-1 ) .* jdlist
   if size(out,2) > 1
      out[:,2:end] .* 2
   end
   return (2s+1)*sum(out)
end
function thermeff(ν,El,kT,Q)
   x = (kT*csl)^-1
   return ν*exp(-El*x)*(1.0-exp(-ν*x)) / Q*csl
end


function jbjklister(jmin,jmax,mΔj)
   out = zeros(0,2)
   for j in jmin:(jmax-1)
      out = vcat(out,[collect(j:(mΔj+j)) fill(j,mΔj+1)])
   end
   out = vcat(out,[collect(jmax:(mΔj+jmax-1)) fill(jmax,mΔj)])
   return out
end
function jbjklisterfull(jmin,jmax,mΔj)
   out = zeros(0,2)
   for j in jmin:jmax
      p1 = collect(max(jmin,(j-mΔj)):(j+mΔj))
      out = vcat(out,[fill(j,length(p1)) p1])
   end
   return out
end
function jbjklister(jmin,jmax,mΔJ,σbσk::Bool)
   if σbσk
      return jbjklister(jmin,jmax,mΔJ)
   else
      return jbjklisterfull(jmin,jmax,mΔJ)
   end
end

function enact_mu(O::Mu,ψb::Psi,σb::Int,ψk::Psi,σk::Int,wvs::Eigs)
   out = O.val
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
         part = eval_tmu(ψb.T, O.tf[i], ψk.T, wvs.top)
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
      out = kron(I(ψk.T.l), out)
   else # !isnothing(wvs.ttp.vals) # two stage no op
      out = kron(wvs.ttp.vecs[:,:,σb]' * wvs.ttp.vecs[:,:,σk],out) # ⟨ ψt' | ψt ⟩
   end #tor if
   # ⟨ ψtr' | O_f | ψtr ⟩
   #@show wvs.rst.vecs[1:size(out,1), jinds(ψb.R), σb]
   #@show wvs.rst.vecs[1:size(out,2), jinds(ψk.R), σk]
   out = #sparsify!( 
      wvs.rst.vecs[1:size(out,1), jinds(ψb.R), σb]' *
      out * 
      wvs.rst.vecs[1:size(out,2), jinds(ψk.R), σk]
   #,1e-5)
   return sparse(out)
end

function μ2_build(Ops::Vector{Mu},ψb::Psi,σb::Int,ψk::Psi,σk::Int,wvs::Eigs)
   μmat = enact_mu(Ops[1], ψb,σb,ψk,σk, wvs)
   for i ∈ 2:length(Ops) # eachindex(Ops)
      μmat += enact_mu(Ops[i], ψb,σb,ψk,σk,wvs)
   end
   μmat .^= 2
   if σb==σk && ψb.R.J==ψk.R.J
      μmat = sparse( UpperTriangular(μmat) - Diagonal(μmat) )
   end

   return μmat
end

function μ_proc!(ctrl,ints,Ops,ψb,σb,ψk,σk,wvs)
   binds = jinds(ψb.R.J, ψb.R.S, ctrl.vtmax+1) 
   kinds = jinds(ψk.R.J, ψk.R.S, ctrl.vtmax+1) 
   rind, cind, tints = findnz( μ2_build(Ops, ψb,σb, ψk,σk, wvs) )
   for l ∈ eachindex(tints)
      ints[binds[rind[l]], kinds[cind[l]] ] = tints[l]
   end
   return sparse(ints)
end


function μσbσk_stage(ctrl,Ops, kT,Q, ϕb,σb, ϕk,σk, wvs)
   mΔj = 1
   ints = spzeros(size(wvs.rst.vals,1), size(wvs.rst.vals,1))
   jbjk = jbjklister(0.5*isodd(ctrl.S), ctrl.Jmax,mΔj, σb==σk)
   for i ∈ 1:size(jbjk,1)
      ψb = Psi(RPsi(jbjk[i,1],ctrl.S), ϕb)
      ψk = Psi(RPsi(jbjk[i,2],ctrl.S), ϕk)
      ints = μ_proc!(ctrl,ints,Ops,ψb,σb,ψk,σk,wvs)
   end
   outfs = zeros(Float64, nnz(ints),3)
   outis = zeros(Int, nnz(ints),4)
   rinds, cinds, vints = findnz(ints)
   for i ∈ 1:nnz(ints)
      r = rinds[i]
      c = cinds[i]
      ν = wvs.rst.vals[r,σb] - wvs.rst.vals[c,σk]
      if ctrl.νrange[2] > ν*1e-3 > ctrl.νrange[1]
         outfs[i,1] = ν
         outfs[i,2] = vints[i]*thermeff(ν,wvs.rst.vals[c,σk], kT,Q)
         outfs[i,3] = wvs.rst.vals[c,σk]/csl
         outis[i,:] = [r σb c σk]
      elseif ctrl.νrange[2] > -ν*1e-3 > ctrl.νrange[1]
         outfs[i,1] = -ν
         outfs[i,2] = vints[i]*thermeff(-ν,wvs.rst.vals[r,σb], kT,Q)
         outfs[i,3] = wvs.rst.vals[r,σb]/csl
         outis[i,:] = [c σk r σb]
      else
      end
   end
   filt = findall(!iszero, outfs[:,1])
   #outis = outis[filt,:] 
   #outfs = outfs[filt,:]
   #filt = findall(x -> abs(x)>ctrl.INTthres, outfs[:,2])
   return outis[filt,:], outfs[filt,:]
end

function tracalc(ctrl,Ops,wvs)
   inds = zeros(Int, 0,4)
   σs = σgen(ctrl.NFOLD)
   frqs = zeros(Float64, 0,3)
   Q = qtot_calc(ctrl.Jmax,ctrl.S,ctrl.TK,wvs)
   kT = ctrl.TK * kb
   for σb ∈ 1:σcount(ctrl.NFOLD)#, σk ∈ 1:σb
      ϕb = TTPsi(ctrl.NFOLD,σs[:,σb],ctrl.mcalc) 
#     ϕk = TTPsi(ctrl.NFOLD,σs[:,σk])
      #ϕk = TTPsi(ctrl.NFOLD,σs[:,σb])
      tinds, tfrqs = μσbσk_stage(ctrl, Ops, kT,Q, ϕb,σb, ϕb,σb, wvs)
      inds = vcat(inds,tinds)
      frqs = vcat(frqs,tfrqs)
   end
   perm = sortperm(frqs[:,1])
   frqs = frqs[perm,:]
   inds = inds[perm,:]
   return frqs, inds
end
