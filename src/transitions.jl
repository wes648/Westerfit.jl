

Tμ(q::Int)::Int = q + 2

function μproc(μf,μo)
   len = maximum(μo)+1
   μo = μo[[3;5;8;9],:]
   csop = sum(μo[3:4,:],dims=1)
   out = zeros(Float64,3,len)
   for i in 1:length(μf)
      #v = μo[:,i]
      #ax = argmax(v[[1;2;4]])
      #ind = csop[i] + 1
      out[argmax(μo[[1;2;4],i]), csop[i] + 1] = μf[i]
   end
   @simd for i in 1:len
      out[:,i] = cart2sphr(out[:,i])
   end
   return out
end

function μred(s::Float64,jb::Float64,nb::Int,jk::Float64,nk::Int)::Float64
   return wig6j(nk,jk,s,
                jb,nb,1)*jnred(jb,nb)*jnred(jk,nk)#*√(2*nb+1)
end
function μelem(pr::Float64,q,s::Float64,jb::Float64,nb,kb,jk::Float64,nk,kk)::Array{Float64,2}
   @. return pr*wig3j( nb,1,nk,
                      -kb,q,kk)*μred(s,jb,nb,jk,nk)*
            powneg1(s+jb-kb+nb+nk)
end
function μmat(μs,s,jb,jk)
   lb = Int((2.0*s+1.0)*(2.0*jb+1.0))
   lk = Int((2.0*s+1.0)*(2.0*jk+1.0))
   nks = permutedims(ngeni(jk,s,lb))
   kks = permutedims(kgeni(jk,s,lb))
   nbs = ngeni(jb,s,lk)
   kbs = kgeni(jb,s,lk)
   out = zeros(lb,lk)
   @simd for q in -1:1
      out += μelem(μs[Tμ(q)],q,s,jb,nbs,kbs,jk,nks,kks)
   end
   return out
end
function μpart(μ,s,jb,nb,kb,jk,nk,kk)
   out = zeros(size(nb))
   @simd for q in -1:1
      out += μelem(μ[Tμ(q)],q,s,jb,nb,kb,jk,nk,kk)
   end
   return out
end
function fullμmat_prev(μs,nf,mcalc,s,jb,σb,jk,σk)
   lb = Int((2.0*s+1.0)*(2.0*jb+1.0))
   lk = Int((2.0*s+1.0)*(2.0*jk+1.0))
   nks = permutedims(ngeni(jk,s,lb))
   kks = permutedims(kgeni(jk,s,lb))
   nbs = ngeni(jb,s,lk)
   kbs = kgeni(jb,s,lk)
   mks = msgen(nf,mcalc,σk)
   mbs = msgen(nf,mcalc,σb)
   mks = kron(mks,ones(Int,1,length(mbs)))
   mbs = permutedims(kron(mbs,ones(Int,1,size(mks,1))))
   out = kron(cosp(0,mbs,mks),μpart(μs[:,1],s,jb,nbs,kbs,jk,nks,kks))
   @simd for i in 2:size(μs,2)
      out += kron(cosp(i-1,mbs,mks),μpart(μs[:,i],s,jb,nbs,kbs,jk,nks,kks))
   end
   return out
end
function fullμmat(μs,nf,mcalc,s,jb,σb,jk,σk)
   lb = Int((2.0*s+1.0)*(2.0*jb+1.0))
   lk = Int((2.0*s+1.0)*(2.0*jk+1.0))
   nks = ngeni(jk,s,lb)
   kks = kgeni(jk,s,lb)
   nbs = permutedims(ngeni(jb,s,lk))
   kbs = permutedims(kgeni(jb,s,lk))
   mks = msgen(nf,mcalc,σk)
   mbs = msgen(nf,mcalc,σb)
   mks = permutedims(kron(mks,ones(Int,1,length(mbs))))
   mbs = kron(mbs,ones(Int,1,size(mks,1)))
   out = kron(cosp(0,mbs,mks),μpart(μs[:,1],s,jb,nbs,kbs,jk,nks,kks))
   @simd for i in 2:size(μs,2)
      out += kron(cosp(i-1,mbs,mks),μpart(μs[:,i],s,jb,nbs,kbs,jk,nks,kks))
   end
   return out
end
function intmat(μ,INTTHRESH,nf,mcalc,s,jb,σb,vecb,jk,σk,veck)
   out = fullμmat(μ,nf,mcalc,s,jk,σk,jb,σb)
   out = sparse(transpose(vecb)*(out)*veck)
   out .*= out
   return droptol!(out,INTTHRESH*0.01)
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

function tracalc_nocat(μ::Array{Float64},kbT,Qrt,ctrl,jmax,
                       rvals,rvecs,rqns,σr,cvals,cvecs,cqns,σc,uncs)
   s = ctrl["S"]
   vtm = ctrl["vtmax"]
   nf = ctrl["NFOLD"]
   mΔj = 1
   #println(μ)
   if σr == σc
      jbjk = jbjklister(0.5*iseven(Int(2*s+1)),jmax,mΔj)
   else
      jbjk = jbjklisterfull(0.5*iseven(Int(2*s+1)),jmax,mΔj)
   end
   rmsd = Int((2*s+1)*(2*ctrl["mcalc"]+1))
   cmsd = Int((2*s+1)*(2*ctrl["mcalc"]+1))
   frqs = spzeros(length(rvals),length(cvals))
   ints = spzeros(length(rvals),length(cvals))
   exprmin = exp(-minimum(rvals)/kbT)
   expcmin = exp(-minimum(cvals)/kbT)
   #elow = spzeros(length(vals),length(vals))
   #build matrix of intensities
   for i in 1:size(jbjk,1)
      jb = jbjk[i,1]
      jk = jbjk[i,2]
      #println("jb=$jb, jk=$jk")
      #needs to be corrected for the vt scheme
      kinds = jvlinds(jk,s,vtm)
      binds = jvlinds(jb,s,vtm)
      #filter vecs
      kvecs = cvecs[1:Int(2*jk+1)*cmsd,kinds]
      bvecs = rvecs[1:Int(2*jb+1)*rmsd,binds]
      #calculate intensities
      μs = intmat(μ,ctrl["INTthres"],nf,ctrl["mcalc"],s,jb,σr,bvecs,jk,σc,kvecs)
      if (jb==jk)&&(σr==σc)
         μs = sparse(UpperTriangular(μs))
      end
      rind, cind, tints = findnz(μs)
      for l in 1:length(tints)
         b = binds[rind[l]]
         k = kinds[cind[l]]
         ints[b,k] = tints[l]
      end
      #ints[binds[cind],kinds[rind]] = tints
   end
   #println(ints)
#I think I can split this function in two here
   #calculate frequencies & filter out ones out of range
   rinds, cinds, = findnz(ints)
   for i in 1:nnz(ints)
      r = rinds[i]
      c = cinds[i]
      ν = rvals[r] - cvals[c]
      frqs[r,c] = ν*(ctrl["νmin"]≤abs(ν*0.001)≤ctrl["νmax"])
   end
   droptol!(frqs,1e-3) #drops all frequencies below 1kHz
   #assign Elow & qns plus sign correction
   rinds, cinds, νs = findnz(frqs)
   len = nnz(frqs) 
   outfrq = zeros(len,4)
   outqns = zeros(Int,len,12)
   for i in 1:len
      ν = νs[i]
      r = rinds[i]
      c = cinds[i]
      if ν > 0.0
         outfrq[i,1] = ν
         outfrq[i,3] = cvals[c] / csl
         outfrq[i,4] = √abs(uncs[r,r]^2 + uncs[c,c]^2 - 2*uncs[r,c])
         thermfact = abs(exp(-cvals[c]/kbT) - exp(-rvals[r]/kbT))/Qrt
         outfrq[i,2] = ints[r,c]*thermfact*(2*s+1)*(σc+1)*100
         outqns[i,1:6] = rqns[r,:]
         outqns[i,7:12] = cqns[c,:]
      elseif ν < 0.0
         outfrq[i,1] = -ν
         outfrq[i,3] = rvals[r] /csl
         outfrq[i,4] = √abs(uncs[r,r]^2 + uncs[c,c]^2 - 2*uncs[r,c])
         thermfact = abs(exp(-rvals[r]/kbT) - exp(-cvals[c]/kbT))/Qrt
         outfrq[i,2] = ints[r,c]*thermfact*(2*s+1)*(σr+1)*100
         outqns[i,1:6] = cqns[c,:]
         outqns[i,7:12] = rqns[r,:]
      end
   end
   #println(outfrq)
   #println(outqns)
   return outfrq, outqns
end

function tracalc_nocat(μ::Nothing,kbT,Qrt,nf,s,jmax,mcalc,vtm,
                       rvals,rvecs,rqns,σr,cvals,cvecs,cqns,σc)
   println("Yikes! Hard to do intensity calculations without dipole components")
   return nothing,nothing
end


function traerrs_bad(J,σu)
   out = zeros(size(J,1))
   @show J
   J = J.^2
   σu = σu.^2
   out = .√sum(J .* σu',dims=2)
   return out
end
function traerrs(J,cov)
   out = J * abs.(cov) * J'
   #return out
   return zeros(size(out))
   #return diag(out)
end

function approxcovar(params,perm)
   out = zeros(length(perm),length(perm))
   for i in 1:length(perm), j in i:length(perm)
      a = perm[i]
      b = perm[i]
      out[i,j] = √(abs(params[a]*params[b]*1e-6))
   end
   return Symmetric(out)
end


function unccalc_no_fit(ctrl,quns,params,scals,stg,ops,pσ,vecs)
   #(ctrl,nlist,ofreqs,uncs,inds,params,scales,cdo,stg,molnam)
   nf = ctrl["NFOLD"]
   s = ctrl["S"]
   vtm = ctrl["vtmax"]
   #inds = qnconv(quns,nf,s,vtm)
   jlist = jlister(inds)
   perm = permdeterm(scals,stg)
   #THIS JACOBIAN IS COMPLETELY WRONG IT IS FOR TRANSITIONS BUT NEED LEVELS
   J = build_jcbn2!([0.0],ops,jlist,inds,ctrl,vecs,params,perm,scals,stg)
   if pσ==nothing
      H = J' * J
      W = I(size(H,1))
      omc = ones(size(J,1))
      pσ = paramunc(H,W,perm,omc)
   end
   uncs = traerrs(J,pσ)
end








