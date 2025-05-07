"""
This is the file that contains the analytic derivative functions
jacobian - 1 stage
jacobian - 2 stage
hessian - jtwj
hessian - briss 1 stage
hessian - briss 2 stage
"""

function sumder(out,j,s,nf,id,prm,stg,ops,ms,qns)
   ind = id+1
   if ind ≤ length(stg)+11
      check = stg[ind-11]
      while check < zero(check)
         pm = prm[ind]
         out .+= enact(ℋ[ind-11],ψ,val,ur,ut)
         ind += 1
         if ind-11 ≤ length(stg)
            check = stg[ind-11]
         else
            check = 0
         end
      end
   end
   return out
end

function derivmat(id,prm,scl,stg,ψ,ℋ)
   if scl[id] ≤ 0 #should this be ≤ 0 ???
   elseif id ≤ 4 #pure rot
      pr = zeros(4)
      pr[id] = 1.0
      out = hrot2(pr,qns)
      out = kron(I(length(ms)),out)
   elseif 5 ≤ id ≤ 8 #spin-rot
      pr = zeros(5)
      pr[id-4] = 1.0
      out = hsr(pr,j,s,qns)
      out = kron(I(length(ms)),out)
   elseif 9 ≤ id ≤ 11 #qua
      pr = zeros(3)
      pr[id-9] = 1.0
      out = hqu(pr,j,s,qns)
      out = kron(I(length(ms)),out)
   else #user def
      out = enact(ℋ[id-11],ψ,val,ur,ut)
      out .= sumder(out,j,s,nf,id,prm,stg,ops,ms,qns)
   end
   return out
end

function anaderiv(prm,scl,stg,id,ops,j,s,nf,ms,qns,vec)
   mat = derivmat(j,s,nf,id,prm,scl,stg,ops,ms,qns)
   out = transpose(vec)*mat*vec
   return diag(out)
end

function derivcalc(jlist,ops,ctrl,perm,vecs,prm,scl,stg)#removed nf call from here, fix in references
   #sd = Int(2*s+1)
   s = ctrl["S"]
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   #jmin = 0.5*iseven(Int(2*s+1))
   #println(jlist)
   #jmax = jlist[end,1]
   #jfd = Int(2*s+1)*Int(sum(2.0 .* collect(Float64,jmin:jmax) .+ 1.0))
   #vtd = Int(ctrl["vtmax"]+1)
   σcnt = σcount(nf)
   derivs = zeros(Float64,size(vecs,2),σcnt,length(perm))
   for sc in 1:σcnt
      #println(sc)
      σ = sc - 1
      msd = Int((2*mcalc+1)*(2s+1))
      #msd = Int(2*s+1)*mcd
      #mstrt, mstop = mslimit(nf,mcalc,σ)
      #jmsd = Int(msd*(2*jmax+1))
      #jsvd = Int(jfd*vtd)
      ms = msgen(nf,mcalc,σ)
      jsublist = jlist[isequal.(jlist[:,2],σ), 1] .* 0.5
      for j in jsublist
         #println(j)
         jd = Int(2.0*j) + 1
         sind, find = jvdest(j,s,ctrl["vtmax"]) 
         qns = qngen(j,s)
         vec = vecs[1:jd*msd,sind:find,sc]
         for i in 1:length(perm)
            pid = perm[i]
            ders = anaderiv(prm,scl,stg,pid,ops,j,s,nf,ms,qns,vec)
            derivs[sind:find,sc,i] = ders#*scl[pid]
         end#perm loop
      end #j loop
   end#σ loop
   return derivs
end#function

function derivcalc_all(ops,ctrl,perm,vecs,prm,scl,stg,σ)
   #all as in all states
   s = ctrl["S"]
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   ms = msgen(nf,mcalc,σ)
   derivs = zeros(Float64,size(vecs,2),length(perm))
   sd = Int(2.0*s+1.0)
   jmin = 0.5*iseven(sd)
   jmax = ctrl["Jmax"]
   jlist = collect(Float64,jmin:jmax)
   msd = Int((2*mcalc+1)*(2s+1))
   @threads for j in jlist
      jd = Int(2.0*j) + 1
      sind, find = jvdest(j,s,ctrl["vtmax"])
      qns = qngen(j,s)
      vec = vecs[1:jd*msd,sind:find]
      for i in 1:length(perm)
         pid = perm[i]
         ders = anaderiv(prm,scl,stg,pid,ops,j,s,nf,ms,qns,vec)
         derivs[sind:find,i] = ders
      end#perm loop
   end#j loop
   return derivs
end#function
