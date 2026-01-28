
function λgen(μ::Float64,er::Float64)::Float64
   ρ = 0.5
   er /= 1000.0 #convert err to GHz
   λ = ρ*μ*er #λF
   λ += (1.0 - ρ)*μ*er/(1+er) #λARC
   return λ
end
function lbmq_gain(β,λ::Float64,jtw,h,omc,nomc)::Float64
   out = 2β' * (λ*Diagonal(h)*β + jtw*omc)
#   out = 2.0*(sum(abs2,omc) - sum(abs2, omc + J*β))
   if out < 0
      println("fucking gain function")
      out = 0.0
   else
#   out = 2.0*(sum(abs2, omc .- nomc)) / out#abs(out)
      out = 2.0*(sum(abs2, omc) - sum(abs2, nomc)) / out#abs(out)
   end
   return out
end
function lbmq_gain2(β,J,omc,nomc)::Float64
   pred = sum(abs2, (J * β) .+ omc)
   actu = sum(abs2,nomc)
   curr = sum(abs2,omc)
   #@show -pred - curr
   #@show actu < curr
   out = (actu - curr) / (-pred - curr)
   #out = (curr - actu) / ( curr - pred)
   return out
end


function jlister(inds)
   #finds all the unique J & σ pairs
   js = vcat(inds[:,1],inds[:,4])
   σs = vcat(inds[:,2],inds[:,5])
   temp = fill((0,0),size(js))
   for i in 1:size(js,1)
      temp[i] = (js[i],σs[i])
   end
   temp = unique(temp)
   jsσs = zeros(Int,size(temp)[1],2)
   jsσs[:,1] = (x->x[1]).(temp)
   jsσs[:,2] = (x->x[2]).(temp)
   jsσs = jsσs[sortperm(jsσs[:,1]),:]
#   jsσs = jsσs[sortperm(jsσs[:,2])]
   return jsσs
end

function rmscalc(vals,inds,ofreqs)
   cfreqs = zero(ofreqs)
   #@show inds
   #@show size(ofreqs)
   #@show size(vals)
   #@show inds[1:5,:]
   #@show inds[end-5:end,:]
#@threads for i in 1:size(cfreqs,1)
   for i in 1:size(cfreqs,1)
      cfreqs[i] = vals[inds[i,3],inds[i,2]+1] - vals[inds[i,6],inds[i,5]+1]
   end
   #println(cfreqs)
   omc = ofreqs - cfreqs
   rms = BLAS.nrm2(omc)/√length(omc)
   #rms = norm(omc)/length(omc)
   return rms, omc, cfreqs
end
function rmscalc(cfreqs,ofreqs)
   omc = ofreqs - cfreqs
   rms = BLAS.nrm2(omc)/√length(omc)
   return rms, omc
end

function sumder(out,j,s,nf,rpid,prm,stg,ops,ms,qns)
   ind = rpid+1
   if ind ≤ length(stg)+18
      check = stg[ind-18]
      while check < zero(check)
         pm = prm[ind]
         out .+= tsr_op(pm,j,s,qns,ms, ops[:, ind-18] )
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

function derivmat(j,s,nf,rpid,prm,scl,stg,ops,ms,qns)
   if scl[rpid] < 0 #should this be ≤ 0 ???
   elseif rpid ≤ 4 #pure rot
      pr = zeros(4)
      pr[rpid] = 1.0
      out = hrot2(pr,qns)
      out = kron(I(length(ms)),out)
   elseif 5 ≤ rpid ≤ 9 #spin-rot
      pr = zeros(5)
      pr[rpid-4] = 1.0
      out = hsr(pr,j,s,qns)
      out = kron(I(length(ms)),out)
   elseif 10 ≤ rpid ≤ 12 #qua
      pr = zeros(3)
      pr[rpid-9] = 1.0
      out = hqu(pr,j,s,qns)
      out = kron(I(length(ms)),out)
   elseif rpid==13 # F
      pr = [1.;0.;0.;0.]
      out = kron(htor2(pr,nf,ms), I(size(qns,1)))
   elseif rpid==16 # Vnf
      pr = [0.;0.;0.;1.]
      out = kron(htor2(pr,nf,ms), I(size(qns,1)))
   elseif rpid==14 # ρzF
      out = kron(pa_op(ms,1), nz_op(qns,1))
   elseif rpid==15 # ρxF
      out = kron(pa_op(ms,1), npm_op(qns,1)) 
   elseif rpid==17 # ηz
      out = kron(pa_op(ms,1), sz_op(j,s,qns,1)) 
   elseif rpid==18 # ηx
      out = kron(pa_op(ms,1), spm_op(j,s,qns,1))
   else #user def
      out = tsr_op(1.0,j,s,qns,ms,ops[:,rpid-18] )
      out .= sumder(out,j,s,nf,rpid,prm,stg,ops,ms,qns)
   end
   return out
end

function anaderiv(prm,scl,stg,rpid,ops,j,s,nf,ms,qns,vec)
   mat = derivmat(j,s,nf,rpid,prm,scl,stg,ops,ms,qns)
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
   #@threads for j in jlist
   for j in jlist
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

function build_jcbn2!(jcbn,ops,jlist,inds,ctrl,vecs,params,perm,scals,stg)
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   #jlist = unique(vcat(inds[:,1:3],inds[:,4:6]))
   jcbn = zeros(Float64,size(inds,1),length(perm))
   deriv = derivcalc(jlist,ops,ctrl,perm,vecs,params,scals,stg)
   #@threads for p in 1:length(perm)
   for p in 1:length(perm)
   for a in 1:size(inds,1)
      jcbn[a,p] = deriv[inds[a,3],inds[a,2]+1,p] - deriv[inds[a,6],inds[a,5]+1,p]
   end
   end
   #@show jcbn
   return jcbn
end
function build_jcbn3!(jcbn,ops,jlist,inds,ctrl,vecs,params,perm,scals,stg)
   #This jacobian builder builds for ALL operators. It is used for predicting the 
   #   new energy levels from the pre-existing eigenvectors
   #nf = ctrl["NFOLD"]
   #println("hi")
   #@show length(params)
   #mcalc = ctrl["mcalc"]
   tperm = collect(1:length(params))[(scals .≥ 0) .* (vcat(ones(18),stg) .> 0)]
   #@show tperm
   #jlist = unique(vcat(inds[:,1:3],inds[:,4:6]))
   jcbn = zeros(Float64,size(inds,1),length(tperm))
   deriv = derivcalc(jlist,ops,ctrl,tperm,vecs,params,scals,stg)
   #@threads 
   for p in 1:length(tperm)
   for a in 1:size(inds,1)
      jcbn[a,p] = deriv[inds[a,3],inds[a,2]+1,p] - deriv[inds[a,6],inds[a,5]+1,p]
   end
   end
   #@show jcbn
   return jcbn
end

function build_jcbn_sim(ops,inds,ctrl,vecs,params,perm,scals,stg)
   #This function is deeply fucked up. The difference shouldn't be there but I'm not sure what
   #  the fucking inds input currently is. I will fix this later
   nf = ctrl["NFOLD"]
   mcalc = ctrl["mcalc"]
   jlist = #copy from tsrcalc
   #jlist = unique(vcat(inds[:,1:3],inds[:,4:6]))
   jcbn = zeros(Float64,size(inds)[1],length(perm))
   deriv = derivcalc(jlist,ops,ctrl,perm,vecs,nf,params,scals,stg)
   #@threads 
   for p in 1:length(perm)
   for a in 1:size(inds,1)
      jcbn[a,p] = deriv[inds[a,3],inds[a,2]+1,p] - deriv[inds[a,6],inds[a,5]+1,p]
   end
   end
   return jcbn
end

function build_hess!(hssn,dk,jtw,jcbn,weights)
   jtw = transpose(jcbn)*weights
   hssn = jtw*jcbn
   #Threads.@threads for i in 1:size(hssn)[1]
   #   dk[i,i] = norm(hssn[:,i])
   #end
   return hssn, dk, jtw
end
function build_hess!(hssn,jtw,jcbn,weights)
   mul!(jtw,jcbn',weights)
   mul!(hssn,jtw,jcbn)
   #return hssn, jtw
end
function build_hess(jtw,jcbn,weights)
   jtw = transpose(jcbn)*weights
   hssn = jtw*jcbn
   return hssn, jtw
end
function hstab(hssn,params)
   for i in 1:length(params)
      check = abs(hssn[i,i] / params[i])
      if check < 1e-12
         temp = abs(params[i])
         hssn[:,i] .*= temp
         hssn[i,:] .*= temp
      end
   end
end

function tsrapprox(j,β,stg)::Vector{Float64}
   out = zeros(size(j))
   out = j* β[collect(1:length(β))[(vcat(ones(18),stg) .> 0)]]
   return out#vec(sum(j,dims=2))
end

function paramunc(H,W,perm,omc)
   uncs = zeros(size(H,1))
   try
      uncs = diag(inv(Symmetric(H)))
   catch 
      uncs = diag(inv(H))
   end
   uncs .*= (omc' * W * omc)/(length(omc)-length(perm))
   return □rt.(uncs)
end
function correl(H)
   out = zeros(size(H))
   for i in 1:size(H,1), j in i:size(H,2)
      out[i,j] = H[i,j] / √(H[i,i]*H[j,j])
   end
   return Symmetric(out)
end
function covarr(corr,pσ)
   out = zeros(size(corr))
   for i in 1:size(corr,1), j in i:size(corr,2)
      out[i,j] = corr[i,j] * pσ[i] * pσ[j]
   end
   return Symmetric(out)
end
function covarr2(hess,omc)
   σ2 = 2*sum(abs2,omc)/(size(omc,1) - size(hess,1))^2
   return σ2 .* inv(hess)
end

function permdeterm(scls,stgs)
   out = collect(1:length(scls))[(scls .> 0) .* (vcat(ones(18),stgs) .> 0)]
end

function outputstart(molnam,λ,rms)
   io = open("$molnam.out", "a")
   println(io,"Initial RMS = $rms MHz")
   println(io,"Initial λ = $λ")
#   println(io,"")
#   println(io,"-------------------------------------")
#   println(io,"")
   close(io)
end


function fincheck!(conv,endp,rms,βf,λlm,goal,check,ϵ0,ϵ1,ϵ2,counter,LIMIT,prms,grad)
   if (rms ≤ goal)#&&(counter > 1)
      println("A miracle has come to pass. The fit has converged")
      endp = "converge"
      conv = true
   elseif (check < ϵ0)
      println("The RMS has stopped decreasing. Hopefully it is low")
      endp = "RMS"
      conv = true
#   elseif (norm(βf))<ϵ1*(norm(prms)+ϵ1)
   elseif norm(βf ./ prms)<ϵ1
   #This stopping criteria needs to be scaled for the wildly varying parameter
   #magnitues. 3 dec 24
      slλ = (@sprintf("%0.4f", log10(λlm)))
      println("It would appear step size has converged. log₁₀(λ) = $slλ")
      @show norm(βf)
      endp = "step size"
      conv = true
   elseif (λlm > 1.0e+9)#&&(Δlm == 0.0)
      println("λlm exceeded threshold.")
      println("If you were using the turducken, try again without it")
      endp = "LMthresh"
      conv = true
   elseif norm(grad) < ϵ2
      println("Gradient is now quite small! This should be good")
      endp = "grad"
      conv = true
   elseif counter ≥ LIMIT
      println("Alas, the iteration count has exceeded the limit")
      endp = "iter"
      conv = true
   else
   end #check if
   return conv, endp
end

function opt_frame(ctrl,nlist,ofreqs,uncs,inds,params,scales,cdo,stg,molnam)
   vals,vecs, = tsrcalc2(params,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
   LIMIT = ctrl["maxiter"]
   paramarray = zeros(Float64, length(params), LIMIT+1)
   paramarray[:,1]=params
   oparams = paramarray[:,1]
   rms, omc, = rmscalc(vals, inds, ofreqs)
   perm = permdeterm(scales,stg)
   println("Initial RMS = $rms")
   goal = BLAS.nrm2(uncs)/√length(uncs)*ctrl["goal"]
   W = Diagonal(1.0 ./ uncs)
   #RHOTHRES = -1.0E-6
   ϵ0 = 0.1E-6
   ϵ1 = 0.1E-6
   μlm = ctrl["λlm0"]#(rms + rms^2)#*0.0
   λlm = λgen(μlm, rms) 
   oλlm = λlm
   println("Initial λ = $λlm")
   Δlm = nrm2(params[perm])/length(perm)
   counter = 0
   BAD = 0

   io = open("$molnam.out", "a")
   println(io,"Initial RMS = $rms MHz")
   println(io,"Initial λ = $λlm")
   #println(io,"")
   #println(io,"-------------------------------------")
   #println(io,"")
   close(io)

   nparams = copy(params)
   βf = zero(perm) #step
   J = zeros(Float64,size(inds,1),length(perm)) #Jacobian
   jtw = zero(omc) #jtwient
   J = build_jcbn2!(J,cdo,nlist,inds,ctrl,vecs,params,perm,scales,stg)
   H, jtw = build_hess(jtw,J,W)
   if true ∈ isnan.(H)
      println("FUCKING FUCKING FUCK. NaN in Hessian")
   end
   endpoint = "not yet"
   converged=false
while converged==false
   λlm = λgen(μlm, rms) 
   β,nomc,nrms,nparams = step_calc()
   check = abs(nrms-rms)/rms
   ρlm = lbmq_gain(βf,λlm,jtw*omc,omc,nomc)
   if nrms < rms
      rms = nrms
      omc = nomc
      params .= nparams
      @time J = build_jcbn2!(J,cdo,nlist,inds,ctrl,vecs,params,perm,scales,stg)
      H, jtw = build_hess(jtw,J,W)
      counter += 1
      paramarray[:,counter+1] = params
      #sρlm = (@sprintf("%0.4f", ρlm))
      srms = (@sprintf("%0.4f", rms))
      slλ = (@sprintf("%0.4f", log10(λlm)))
      #sΔ = (@sprintf("%0.6f", Δlm))
      scounter = lpad(counter,3)
      println("After $scounter iterations, RMS = $srms, log₁₀(λ) = $slλ")
      iterationwriter(molnam,paramarray,srms,scounter,slλ,βf,perm)
      μlm /= 30.0
   else
      μlm = max(4.0*μlm,1.0E-24)
   end
   if (rms ≤ goal)#&&(counter > 1)
      println("A miracle has come to pass. The fit has converged")
      endpoint = "converge"
      break
   elseif (check < ϵ0)
      println("The RMS has stopped decreasing. Hopefully it is low")
      endpoint = "RMS"
      break
   elseif (norm(βf))<ϵ1*(norm(params[perm])+ϵ1)
      slλ = (@sprintf("%0.4f", log10(λlm)))
      println("It would appear step size has converged. log₁₀(λ) = $slλ")
      endpoint = "step size"
      break
   elseif (λlm > 1.0e+9)&&(Δlm == 0.0)
      println("λlm exceeded threshold.")
      println("If you were using the turducken, try again without it")
      endpoint = "LMthresh"
      break
   elseif counter ≥ LIMIT
      println("Alas, the iteration count has exceeded the limit")
      endpoint = "iter"
      break
   else
   end #check if
end#while
   frms, fomc, fcfrqs = rmscalc(vals, inds, ofreqs)
   puncs = zeros(size(params))
   puncs[perm] = paramunc(H,W,perm,omc)
   covarmat = covarr(correl(H),puncs)
   params[1:18], puncs[1:18] = fullrecov(params[1:18],puncs[1:18],ctrl["Irrep"])
   slλ = (@sprintf("%0.4f", log10(λlm)))
   outputfinal(molnam,ctrl,frms,counter,slλ,puncs,params,endpoint)
   if ctrl["overwrite"]==true
      println("Writing new input file at $molnam.inp. Previous file has moved to $molnam","1.inp")
      inpwriter(molnam, params, scales)
   end
   return params, covarmat, fomc, fcfrqs, vals
end
