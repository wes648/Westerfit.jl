
Tμ(q::Int)::Int = q + 2

function μproc(μf,μo)
   len = maximum(μo)+1
   axop = μo[1:3,:]
   csop = sum(μo[7:8,:],dims=1)
   out = zeros(Float64,3,len)
   for i in 1:length(μf)
      ax = argmax(axop[:,i])
      ind = csop[i] + 1
      out[ax,ind] = μf[i]
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
               (-1)^(s+jb-kb+nb+nk)
end
function μmat(μs,s,jb,jk)
   lb = Int((2.0*s+1.0)*(2.0*jb+1.0))
   lk = Int((2.0*s+1.0)*(2.0*jk+1.0))
   nks = Matrix(transpose(ngeni(jk,s,lb)))
   kks = Matrix(transpose(kgeni(jk,s,lb)))
   nbs = ngeni(jb,s,lk)
   kbs = kgeni(jb,s,lk)
   out = zeros(lb,lk)
   @simd for q in -1:1
      out += μelem(μs[T(1,q)],q,s,jb,nbs,kbs,jk,nks,kks)
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
   nks = Matrix(transpose(ngeni(jk,s,lb)))
   kks = Matrix(transpose(kgeni(jk,s,lb)))
   nbs = ngeni(jb,s,lk)
   kbs = kgeni(jb,s,lk)
   mks = msbuilder(nf,mcalc,σk)
   mbs = msbuilder(nf,mcalc,σb)
   mks = kron(mks,ones(Int,1,length(mbs)))
   mbs = Matrix(transpose(kron(mbs,ones(Int,1,size(mks,1)))))
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
   nbs = Matrix(transpose(ngeni(jb,s,lk)))
   kbs = Matrix(transpose(kgeni(jb,s,lk)))
   mks = msbuilder(nf,mcalc,σk)
   mbs = msbuilder(nf,mcalc,σb)
   mks = Matrix(transpose(kron(mks,ones(Int,1,length(mbs)))))
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
                       rvals,rvecs,rqns,σr,cvals,cvecs,cqns,σc)
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
   rmsd = Int((2*s+1)*(2*ctrl["mcalc"]+1+(σtype(nf,σr)==2)))
   cmsd = Int((2*s+1)*(2*ctrl["mcalc"]+1+(σtype(nf,σc)==2)))
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
   #calculate frequencies & filter out ones out of range
   rinds, cinds, = findnz(ints)
   for i in 1:nnz(ints)
      r = rinds[i]
      c = cinds[i]
      ν = rvals[r] - cvals[c]
      frqs[r,c] = ν*(ctrl["νmin"]≤abs(ν*0.001)≤ctrl["νmax"])
   end
   dropzeros!(frqs)
   #assign Elow & qns plus sign correction
   rinds, cinds, νs = findnz(frqs)
   len = nnz(frqs)
   outfrq = zeros(len,3)
   outqns = zeros(Int,len,12)
   for i in 1:len
      ν = νs[i]
      r = rinds[i]
      c = cinds[i]
      if ν > 0.0
         outfrq[i,1] = ν
         outfrq[i,3] = cvals[c] / csl
         thermfact = abs(-exp(-cvals[c]/kbT) + exprmin)/Qrt
         outfrq[i,2] = ints[r,c]#*thermfact
         outqns[i,1:6] = rqns[r,:]
         outqns[i,7:12] = cqns[c,:]
      elseif ν < 0.0
         outfrq[i,1] = -ν
         outfrq[i,3] = rvals[r] /csl
         thermfact = abs(-exp(-rvals[r]/kbT) + expcmin)/Qrt
         outfrq[i,2] = ints[r,c]#*thermfact
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
   println("Yikes! Hard to do intensity calculations with out dipole components")
   return nothing,nothing
end
