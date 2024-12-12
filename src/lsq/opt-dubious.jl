"""
Apparently I have to make a new fucking optimizer because some asshole decided 
   to make a new version of westerfit that uses a two stage diagonalization
   process and the data structures won't line up
functions needed:
1. new derivative element calculator -- done?
2. new jacobian calculator -- done?
3. new lbmq implementation to support new data structures
"""

function sumder_2stg(out,j,s,nf,rpid,prm,stg,ops,ms,qns,tvecs)
   ind = rpid+1
   if ind ≤ length(stg)+18
      check = stg[ind-18]
      while check < zero(check)
         pm = prm[ind]
         out .+= twostg_op(pm,j,s,qns,ms, ops[:, ind-18], tvecs )
         #if sum(ms .% 3)≈0 && j==1.5
         #   @show out
         #end
         ind += 1
         if ind-18 ≤ length(stg)
            check = stg[ind-18]
         else
            check = 0
         end
      end
   end
   return out
end
function derivmat_2stg(j,s,nf,rpid,prm,scl,stg,ops,ms,qns,tvecs,mmax)
   if scl[rpid] < 0 #should this be ≤ 0 ???
   elseif rpid ≤ 4 #pure rot
      pr = zeros(4)
      pr[rpid] = 1.0
      out = hrot2(pr,qns)
      out = kron(I(mmax+1),out)
   elseif 5 ≤ rpid ≤ 9 #spin-rot
      pr = zeros(5)
      pr[rpid-4] = 1.0
      out = hsr(pr,j,s,qns)
      out = kron(I(mmax+1),out)
   elseif 10 ≤ rpid ≤ 12 #qua
      pr = zeros(3)
      pr[rpid-9] = 1.0
      out = hqu(pr,j,s,qns)
      out = kron(I(mmax+1),out)
   elseif rpid==13 # F
      pr = [1.;0.;0.;0.]
      out = kron(tvecs' * htor2(pr,nf,ms) * tvecs,  I(size(qns,1)))
   elseif rpid==16 # Vnf
      pr = [0.;0.;0.;1.]
      out = kron(tvecs' * htor2(pr,nf,ms) * tvecs,  I(size(qns,1)))
   elseif rpid==14 # ρzF
      out = kron(tvecs' * pa_op(ms,1) * tvecs, nz_op(qns,1))
   elseif rpid==15 # ρxF
      out = kron(tvecs' * pa_op(ms,1) * tvecs, npm_op(qns,1)) 
   elseif rpid==17 # ηz
      out = kron(tvecs' * pa_op(ms,1) * tvecs, sz_op(j,s,qns,1)) 
   elseif rpid==18 # ηx
      out = kron(tvecs' * pa_op(ms,1) * tvec, spm_op(j,s,qns,1))
   else #user def
      out = twostg_op(1.0,j,s,qns,ms,ops[:,rpid-18], tvecs )
      out .= sumder_2stg(out,j,s,nf,rpid,prm,stg,ops,ms,qns,tvecs)
   end
   return out
end
function anaderiv_2stg(prm,scl,stg,rpid,ops,j,s,nf,ms,qns,vec,tvec,mmax)
   mat = derivmat_2stg(j,s,nf,rpid,prm,scl,stg,ops,ms,qns,tvec,mmax)
   out = transpose(vec)*mat*vec
   return diag(out)
end
function derivcalc_2stg(jlist,ops,ctrl,perm,vecs,prm,scl,stg,tvecs)
   s = ctrl["S"]
   nf = ctrl["NFOLD"]
   σcnt = σcount(nf)
   derivs = zeros(Float64,size(vecs,2),σcnt,length(perm))
   msd = Int((ctrl["mmax"]+1)*(2s+1))
   for sc in 1:σcnt
      σ = sc - 1
      ms = msgen(nf,ctrl["mcalc"],σ)
      tvcs = tvecs[:,:,sc]
      jsublist = jlist[isequal.(jlist[:,2],σ), 1] .* 0.5
      for j in jsublist
         #println(j)
         jd = Int(2.0*j) + 1
         sind, find = jvdest(j,s,ctrl["vtmax"]) 
         qns = qngen(j,s)
         vec = vecs[1:jd*msd,sind:find,sc]
         for i in 1:length(perm)
            pid = perm[i]
            ders = anaderiv_2stg(prm,scl,stg,pid,ops,j,s,nf,ms,qns,vec,tvcs,ctrl["mmax"])
            derivs[sind:find,sc,i] = ders#*scl[pid]
         end#perm loop
      end #j loop
   end#σ loop
   return derivs
end#function

#function derivmat_2stg(j,s,nf,rpid,ms,qns,tvecs,mmax)
#   if rpid ≤ 4 #pure rot
#      pr = zeros(4)
#      pr[rpid] = 1.0
#      out = hrot2(pr,qns)
#      out = kron(I(mmax+1),out)
#   elseif 5 ≤ rpid ≤ 9 #spin-rot
#      pr = zeros(5)
#      pr[rpid-4] = 1.0
#      out = hsr(pr,j,s,qns)
#      out = kron(I(mmax+1),out)
#   elseif 10 ≤ rpid ≤ 12 #qua
#      pr = zeros(3)
#      pr[rpid-9] = 1.0
#      out = hqu(pr,j,s,qns)
#      out = kron(I(mmax+1),out)
#   elseif rpid==13 # F
#      pr = [1.;0.;0.;0.]
#      out = kron(tvecs' * htor2(pr,ms) * tvecs,  I(size(qns,1)))
#   elseif rpid==16 # Vnf
#      pr = [0.;0.;0.;1.]
#      out = kron(tvecs' * htor2(pr,ms) * tvecs,  I(size(qns,1)))
#   elseif rpid==14 # ρzF
#      out = kron(tvecs' * pa_op(ms,1) * tvecs, nz_op(qns,1))
#   elseif rpid==15 # ρxF
#      out = kron(tvecs' * pa_op(ms,1) * tvecs, npm_op(qns,1)) 
#   elseif rpid==17 # ηz
#      out = kron(tvecs' * pa_op(ms,1) * tvecs, sz_op(j,s,qns,1)) 
#   elseif rpid==18 # ηx
#      out = kron(tvecs' * pa_op(ms,1) * tvec, spm_op(j,s,qns,1))
#   else
#      out = Diagonal(I(size(qns,1)))
#   end
#   return out
#end
function anaderiv_2stg(rpid,j,s,nf,ms,qns,vec,tvec,mmax)
   mat = derivmat_2stg(j,s,nf,rpid,ms,qns,tvec,mmax)
   out = transpose(vec)*mat*vec
   return diag(out)
end
function derivcalc_2stg(jlist,ctrl,perm,vecs,tvecs)
   s = ctrl["S"]
   nf = ctrl["NFOLD"]
   σcnt = σcount(nf)
   derivs = zeros(Float64,size(vecs,2),σcnt,length(perm))
   msd = Int((ctrl["mmax"]+1)*(2s+1))
   for sc in 1:σcnt
      σ = sc - 1
      ms = msgen(nf,ctrl["mcalc"],σ)
      tvcs = tvecs[:,:,sc]
      jsublist = jlist[isequal.(jlist[:,2],σ), 1] .* 0.5
      for j in jsublist
         #println(j)
         jd = Int(2.0*j) + 1
         sind, find = jvdest(j,s,ctrl["vtmax"]) 
         qns = qngen(j,s)
         vec = vecs[1:jd*msd,sind:find,sc]
         for i in 1:length(perm)
            pid = perm[i]
            ders = anaderiv_2stg(pid,j,s,nf,ms,qns,vec,tvcs,ctrl["mmax"])
            derivs[sind:find,sc,i] = ders#*scl[pid]
         end#perm loop
      end #j loop
   end#σ loop
   return derivs
end#function

function derivcalc_2stg_all(ops,ctrl,perm,vecs,prm,scl,stg,σ)
   #all as in all states
   s = ctrl["S"]
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   ms = msgen(nf,mcalc,σ)
   tvcs = tvecs[:,:,σ+1]
   derivs = zeros(Float64,size(vecs,2),length(perm))
   sd = Int(2.0*s+1.0)
   jmin = 0.5*iseven(sd)
   jmax = ctrl["Jmax"]
   jlist = collect(Float64,jmin:jmax)
   msd = Int((ctrl["mmax"]+1)*(2s+1))
   @threads for j in jlist
      jd = Int(2.0*j) + 1
      #tvcs = tvecs[:,:]
      sind, find = jvdest(j,s,ctrl["vtmax"])
      qns = qngen(j,s)
      vec = vecs[1:jd*msd,sind:find]
      for i in 1:length(perm)
         pid = perm[i]
         ders = anaderiv_2stg(prm,scl,stg,pid,ops,j,s,nf,ms,qns,vec,tvcs)
         derivs[sind:find,i] = ders
      end#perm loop
   end#j loop
   return derivs
end#function
function build_jcbn_2stg!(jcbn,ops,jlist,inds,ctrl,vecs,params,perm,scals,stg,tvecs)
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   #jcbn = zeros(Float64,size(inds,1),length(perm))
   deriv = derivcalc_2stg(jlist,ops,ctrl,perm,vecs,params,scals,stg,tvecs)
   td = round.(deriv,sigdigits=5)
   @threads for p in 1:length(perm)
   @simd for a in 1:size(inds,1)
      jcbn[a,p] = deriv[inds[a,3],inds[a,2]+1,p] - deriv[inds[a,6],inds[a,5]+1,p]
   end
   end
   #@show jcbn
   return jcbn
end
function out_jcbn_2stg(molnam::String,ops,jlist,inds,ctrl,vecs,params,scals,stg,tvecs)
   perm = collect(1:length(params))[(params .> 0) .* (vcat(ones(18),stgs) .> 0)]
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   jcbn = zeros(Float64,size(inds,1),length(perm))
   deriv = derivcalc_2stg(jlist,ops,ctrl,perm,vecs,params,scals,stg,tvecs)
   td = round.(deriv,sigdigits=5)
   @threads for p in 1:length(perm)
   @simd for a in 1:size(inds,1)
      jcbn[a,p] = deriv[inds[a,3],inds[a,2]+1,p] - deriv[inds[a,6],inds[a,5]+1,p]
   end
   end
   writedlm("$molnam.der",jcbn,';')
   #@show jcbn
   return jcbn
end
function lbmq_2stp!(H,jtw,omc,λ,β,prms)
   #A = Hermitian(H + λ*Diagonal(abs.(prms) .* diag(H)) )
   A = Hermitian(H + λ*Diagonal(H))
   #A = Hermitian(H + λ*I )
   #A = Hermitian(H + λ*Diagonal(abs.(prms)))
   while isposdef(A)==false #this could be tidier
      λ = max(2.0*λ,1.0E-24)
      A = Hermitian(H + λ*Diagonal(H))
      if isinf(λ); break; end
   end
   if isinf(λ)
      @warn "LB-MQ Matrix not pos-def!"
      println("Make sure you aren't trying to optimize a parameter with value of 0.0.")
      println("This code is about to crash")
   end
   A = cholesky!(A)
   β .= ldiv!(β, A, jtw*omc)# .* scls[perm]
   return β,λ
end
function dogleg_corr!(βlm,Δ,jtw,omc)
   #This dogleg function will only run if βlm is outside trust region
   βsd = -jtw*omc
   if norm(βsd) > Δ #both outside trust region
      βlm .= Δ/norm(βsd) * βsd
   else #only GN out of trust regeion
      #resolve t's simplification for weighted
      t = norm(βsd' * jtw * omc) / norm(J*βsd)^2
      s = (Δ - t*norm(βsd))/norm(βlm - βsd)
      #figure out cleaner version to adjust s so below is less than Δ
      #temp = t*βsd + s*(βlm - t*βsd)
      temp = s*βlm + t*(1-s)*βsd
      while norm(temp) > Δ
         s *= .8
         temp = s*βlm + t*(1-s)*βsd
      end#while
      βlm .= temp
   end#end
   return βlm
end


function lbmq_2stg(ctrl,nlist,ofreqs,uncs,inds,params,scales,cdo,stg,molnam)
   STAGE = ctrl["stages"]
   if STAGE==1
   vals,vecs, = tsrcalc2(params,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
   elseif STAGE==2
   vals,vecs,tvcs, = twostg_calc2(params,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
   else
   println("Someday....")
   end
   GEO = true
   BOLD = 1
   LIMIT = ctrl["maxiter"]

   paramarray = zeros(Float64, length(params), LIMIT+1)
   paramarray[:,1]=params
   oparams = paramarray[:,1]

   rms, omc, = rmscalc(vals, inds, ofreqs)
   lrms = rms
   perm = permdeterm(scales,stg)
   W = Diagonal(1.0 ./ uncs)^2
   wrms = √(omc' *W*omc ./ length(omc))
   goal = BLAS.nrm2(uncs)/√length(uncs)*ctrl["goal"]
   println("Initial  RMS = $rms")
   println("Initial WRMS = $wrms")
   ϵ0 = 0.1E-5 #rms change threshold
   ϵ1 = 0.1E-4 #step size threshold
   ϵ2 = 0.1E-3 #gradient threshold
   μlm = ctrl["λlm0"]#(rms + rms^2)#*0.0
   λlm = λgen(μlm, rms) 
   Δlm = 1.0
   println("Initial λ = $λlm")
   counter = 0
   θ = 0.0
   bad = 0
   outputstart(molnam,λlm,rms)

   nparams = deepcopy(params)
   #puncs = zero(perm)
   βf = zeros(length(perm),2) #step
   βo = zeros(length(perm))
   J = zeros(Float64,size(inds,1),length(perm)) #Jacobian
   jtw = zeros(Float64,length(perm),size(inds,1))
   H = zeros(length(perm),length(perm))
   if STAGE==1
   build_jcbn2!(J,cdo,nlist,inds,ctrl,vecs,params,perm,scales,stg)
   elseif STAGE==2
   build_jcbn_2stg!(J,cdo,nlist,inds,ctrl,vecs,params,perm,scales,stg,tvcs)
   else
   println("Someday....")
   end
   #J, w, omc = linereject(J,W,omc,uncs,ctrl["REJECT"])
   build_hess!(H,jtw,J,W)
   if true ∈ isnan.(H)
      println("FUCKING FUCKING FUCK. NaN in Hessian")
   end
   #uncs = paramunc(uncs,H,perm,omc)
   endpoint = "not yet"
   converged=false
   while converged==false
      λlm = λgen(μlm, rms)
      βf[:,1],λlm = lbmq_2stp!(H,jtw,omc,λlm,βf[:,1],params[perm])
      if GEO  && wrms < 2e2 #only uses accelerator when wrms isn't outragous
         dderiv = sum(J*βf*βf'*jtw,dims=2)[:]
         #dderiv = diag(J*βf*βf'*jtw)
         βf[:,2],λlm = lbmq_2stp!(H,jtw,dderiv,λlm,βf[:,2],params[perm])
         βf[:,2] .*= -0.5
         if abs(norm(βf[:,2])/norm(βf[:,1])) > 0.5 #This is half the threshold factor
            #they suggest 0.75 (so 0.375 here)
            βf[:,2] .*= 0.25*abs(norm(βf[:,1])/norm(βf[:,2]))
         end
      end
      #βf .*= scales[perm]
      β = sum(βf,dims=2)[:]
      β .*= scales[perm]
      @show β
      #if norm(H^(-.5)*β) > Δlm
      if norm(β ./params[perm]) > Δlm
      #   @show H^(-.5)*β
      #   β = dogleg_corr!(β,Δlm,jtw,omc)
         β *= 0.5*Δlm/norm(β ./ params[perm])
      end

      if counter > 0
         θ = dot(β,βo) / (norm(β)*norm(βo))
         θ = round(θ; digits=3)
         @show θ
      end

      nparams[perm] .+= β
      if STAGE==1
      vals,vecs, = tsrcalc2(params,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
      elseif STAGE==2
      nvals,nvecs,tvcs, = twostg_calc2(nparams,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
      else
         println("Someday....")
      end
      nrms, nomc, = rmscalc(nvals,inds,ofreqs)
      nwrms = √(nomc' *W*nomc ./ length(nomc))
      ρlm = lbmq_gain(β,λlm,jtw,H,omc,nomc)
      check = abs(nrms-rms)/rms
      println("ρlm = $ρlm, nrms = $nrms")
      println("Δlm = $Δlm, wrms = $nwrms")
      println("λlm = $λlm, norm(β) = $(norm(β))")
      println("norm(β ./ params) = $(norm(β ./ params[perm]))")
      #if (nrms<rms)&&(ρlm>1e-3)
      #stepcheck = (ρlm>1e-5)
      if (BOLD == 0)
         stepcheck = ((nrms<rms)&&(ρlm>1e-6))
      else
      stepcheck = ((nrms<rms)&&(ρlm>1e-6)) || ((nrms*(1-θ)^BOLD)<0.2*lrms)
      end

      #if ((nrms*(1-θ)^BOLD)< lrms)&&BOLD>0|| ((nrms<rms)&&(ρlm>1e-3))# || bad > 2
      if stepcheck
      #@show θ
      #@show round.(bf, sigdigits=5)
      #@show round.(bf./nparams[perm], sigdigits=5)
         counter += 1
         bad = 0
         rms = nrms
         lrms = min(lrms,rms)
         omc = nomc
         tΔ = round.((nvals .- vals)./sum(βf), sigdigits=5)
         vals .= nvals
         params .= nparams
         vecs .= nvecs
         wrms = √(omc' *W*omc ./ length(omc))
         #println(params)
         #@time J = build_jcbn!(J,cdo,inds,S,ctrl,vecs,params,perm,scales)
      if STAGE==1
      @time J = build_jcbn2!(J,cdo,nlist,inds,ctrl,nvecs,nparams,perm,scales,stg)
      elseif STAGE==2
      @time J = build_jcbn_2stg!(J,cdo,nlist,inds,ctrl,nvecs,nparams,perm,scales,stg,tvcs)
      else
         println("Someday....")
      end      
         #J, w = linereject(J,W,omc,uncs,ctrl["REJECT"])
         #build_hess!(H,jtw,J,W)
         H, jtw = build_hess(jtw,J,W)
         βo = β
            #println(diag(H))
         paramarray[:,counter+1] = params
         #sρlm = (@sprintf("%0.4f", ρlm))
         srms = (@sprintf("%0.4f", rms))
         slλ = (@sprintf("%0.4f", log10(λlm)))
         scounter = lpad(counter,3)
         println("After $scounter iterations, RMS = $srms, log₁₀(λ) = $slλ")
         iterationwriter(molnam,paramarray,rms,counter,λlm,sum(βf,dims=2),perm)
         @show wrms
         #μlm /= 20.0
         λlm /= 20.0
         Δlm *= min(1.1,5.0)
         GEO=true
      else
         #GEO = false
         GEO = GEO != true #inverts GEO
         μlm = max(10.0*μlm,1.0E-24)
         Δlm = max(Δlm*.8,1.0E-24)
         #λlm = max(10.0*λlm,1.0E-24)
         bad += 1
      end
   converged,endpoint = fincheck!(converged,endpoint,rms,βf,
         λlm,goal,check,ϵ0,ϵ1,ϵ2,counter,LIMIT,params[perm],-2jtw*omc)
   end#while
   frms, fomc, fcfrqs = rmscalc(vals, inds, ofreqs)
   puncs = zeros(size(params))
   puncs[perm] = paramunc(H,W,perm,omc)
   covarmat = covarr(correl(H),puncs)
   #covarmat = inv(H)
   #params[1:15] .= paramrecov(params[1:15])
   #uncs[1:15] .= uncrecov(uncs[1:15],params[1:15])
   params[1:18], puncs[1:18] = fullrecov(params[1:18],puncs[1:18],ctrl["Irrep"])
   slλ = (@sprintf("%0.4f", log10(λlm)))
   outputfinal(molnam,ctrl,frms,counter,slλ,puncs,params,endpoint)
   if occursin("J", ctrl["RUNmode"])&&STAGE==2
      out_jcbn_2stg(molnam,ops,nlist,inds,ctrl,vecs,params,scales,stg,tvcs)
   end
   if ctrl["overwrite"]==true
      println("Writing new input file at $molnam.inp. Previous file has moved to $molnam","1.inp")
      inpwriter(molnam, params, scales)
   end
   return params, covarmat, fomc, fcfrqs, vals
end
