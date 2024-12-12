
"""
One day I'll need to rewrite both intensity calculators
"""

function fullμmat_2stg(μs,nf,mcalc,s,jb,σb,jk,σk,tvecs)
   tvk = tvecs[:,:,σk+1]
   tvb = tvecs[:,:,σb+1]
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
   out = kron(tvb' * cosp(0,mbs,mks) * tvk,
              μpart(μs[:,1],s,jb,nbs,kbs,jk,nks,kks))
   @simd for i in 2:size(μs,2)
      out += kron(tvb' * cosp(i-1,mbs,mks) * tvk,
                  μpart(μs[:,i],s,jb,nbs,kbs,jk,nks,kks)) 
   end
   return out
end
function intmat2stg(μ,INTTHRESH,nf,mcalc,s,jb,σb,vecb,jk,σk,veck,tvecs)
   out = fullμmat_2stg(μ,nf,mcalc,s,jk,σk,jb,σb,tvecs)
   out = sparse(transpose(vecb)*(out)*veck)
   out .*= out
   return droptol!(out,INTTHRESH*0.01)
end

function tracalc_twostg(μ::Array{Float64},kbT,Qrt,ctrl,jmax,
                       rvals,rvecs,rqns,σr,cvals,cvecs,cqns,σc,uncs,tvecs)
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
   rmsd = Int((2*s+1)*(ctrl["mmax"]+1))
   cmsd = Int((2*s+1)*(ctrl["mmax"]+1))
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
μs = intmat2stg(μ,ctrl["INTthres"],nf,ctrl["mcalc"],s,jb,σr,bvecs,jk,σc,kvecs,tvecs)
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
         #outfrq[i,4] = √abs(uncs[r,r]^2 + uncs[c,c]^2 - 2*uncs[r,c])
         thermfact = abs(exp(-cvals[c]/kbT) - exp(-rvals[r]/kbT))/Qrt
         outfrq[i,2] = ints[r,c]*thermfact*(2*s+1)*(σc+1)*100
         outqns[i,1:6] = rqns[r,:]
         outqns[i,7:12] = cqns[c,:]
      elseif ν < 0.0
         outfrq[i,1] = -ν
         outfrq[i,3] = rvals[r] /csl
         #outfrq[i,4] = √abs(uncs[r,r]^2 + uncs[c,c]^2 - 2*uncs[r,c])
         thermfact = abs(exp(-rvals[r]/kbT) - exp(-cvals[c]/kbT))/Qrt
         outfrq[i,2] = ints[r,c]*thermfact*(2*s+1)*(σr+1)*100
         outqns[i,1:6] = cqns[c,:]
         outqns[i,7:12] = rqns[r,:]
      end
   end
   filt = outfrq[:,2] .≥ ctrl["INTthres"]
   outfrq = outfrq[filt,:]
   outqns = outqns[filt,:]
   #println(outfrq)
   #println(outqns)
   return outfrq, outqns
end
