"""
Apparently I have to make a new fucking optimizer because some asshole decided 
   to make a new version of westerfit that uses a two stage diagonalization
   process and the data structures won't line up
functions needed:
1. new derivative element calculator -- done?
2. new jacobian calculator -- done?
3. new lbmq implementation to support new data structures
"""

function lbmq_2stp!(H,jtw,omc,λ,β,prms)
   #@show H
   A = Hermitian(H + λ*Diagonal(H))
   while isposdef(A)==false #this could be tidier
      λ = max(2.0*λ,1.0E-24)
      A = Hermitian(H + λ*Diagonal(H))
      if isinf(λ); break; end
   end
   if isinf(λ)
      @warn "LB-MQ Matrix not pos-def!
      Make sure you aren't trying to optimize a parameter with value of 0.0.
      This code is about to crash"
   end
   A = cholesky!(A)
   β .= ldiv!(β, A, -jtw*omc)# .* scls[perm]
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


function lbmq(ctrl,nlist,ofreqs,uncs,inds,params,scales,ℋ,stg,molnam)
   STAGE = ctrl.stages
   @show params
   σs = σgen_indef(ctrl.NFOLD)
   σcnt = maximum(size(σs))
   mcd = Int(2*ctrl.mcalc+1)

#I need to make the sublist

   if ctrl.stages==1
      sd = Int(2*ctrl.S+1)
      vtd = ctrl.vtmax + 1
      jfd = sd*Int(sum(2.0 .* unique(nlist[:,1]) .+ 1.0))
      vals = zeros(Float64,jfd*vtd,σcnt)
      vecs = zeros(Float64,Int(sd*(2*maximum(nlist[:,1])+1)*mcd),jfd*vtd,σcnt)
      @time vals,vecs = tsrcalc_1stg!(vals,vecs,nlist,σs,ctrl,params,stg,ℋ)
   elseif ctrl.stages==2
      sd = Int(2*ctrl.S+1)
      vtd = ctrl.vtmax + 1
      jfd = sd*Int(sum(2.0 .* unique(nlist[:,1]) .+ 1.0))
      vals = zeros(Float64,jfd*vtd,σcnt)
      vl = sd*(2*maximum(nlist[:,1])+1)*(ctrl.mmax+1)
      vecs = zeros(Float64,Int(vl),jfd*vtd,σcnt)
      #initialize tvecs
      tvecs = zeros(2*ctrl.mcalc+1,ctrl.mmax+1,σcnt)
@time vals,vecs,tvals,tvecs = tsrcalc_2stg!(vals,vecs,tvals,tvecs,nlist,σs,ctrl,params,stg,ℋ)
   else
      @warn "Invalid stages number"
   end
   GEO = false#true
   BOLD = 1
   LIMIT = ctrl.maxiter

   nvals = zero(vals)
   nvecs = zero(vecs)

   paramarray = zeros(Float64, length(params), LIMIT+1)
   paramarray[:,1]=params
   oparams = paramarray[:,1]

   rms, omc, = rmscalc(vals, inds, ofreqs)
   lrms = rms
   perm = permdeterm(scales,stg)
   W = Diagonal(1.0 ./ uncs)^2
   #W ./= maximum(W)
   wrms = √(omc' *W*omc ./ length(omc))
   goal = BLAS.nrm2(uncs)/√length(uncs)*ctrl.goal
   println("Initial  RMS = $rms")
   println("Initial wRMS = $wrms")
   ϵ0 = 0.1E-8 #rms change threshold
   ϵ1 = 0.1E-6 #step size threshold
   ϵ2 = 0.1E-3 #gradient threshold
   μlm = ctrl.λlm0#(rms + rms^2)#*0.0
   λlm = λgen(μlm, rms) 
   Δlm = 0.5
   println("Initial λ = $λlm")
   counter = 0
   θ = 0.0
   outputstart(molnam,λlm,rms)

   nparams = deepcopy(params)

   βf = zeros(length(perm),2) #step
   βo = zeros(length(perm))
   J = zeros(Float64,size(inds,1),length(perm)) #Jacobian
   jtw = zeros(Float64,length(perm),size(inds,1))
   H = zeros(length(perm),length(perm))
   if STAGE==1
      println("welcome to 1 stage")
   #J = build_jcbn!(J,cdo,nlist,inds,ctrl,vecs,params,perm,scales,stg)
      build_jcbn!(J,inds,nlist,ℋ,ctrl,perm,vecs,params,scales,stg)
   elseif STAGE==2
      build_jcbn_2stg!(J,cdo,nlist,inds,ctrl,vecs,params,perm,scales,stg,tvcs)
   else
      println("Someday....")
   end
   #jtw = J' * W
   #J, w, omc = linereject(J,W,omc,uncs,ctrl.REJECT)
   build_hess!(H,jtw,J,W)
   if true ∈ isnan.(H)
      println("FUCKING FUCKING FUCK. NaN in Hessian")
   end
   @show H
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
         end#if
      end#GEO
      βf .*= scales[perm]
      β = 10 .*sum(βf,dims=2)[:]
      β .*= scales[perm]
      #if norm(H^(-.5)*β) > Δlm
      if norm(β ./params[perm]) > Δlm
      ##   @show H^(-.5)*β
      ##   β = dogleg_corr!(β,Δlm,jtw,omc)
         β .*= 0.5*Δlm/norm(β ./ params[perm])
      end

      if counter > 0
         θ = dot(β,βo) / (norm(β)*norm(βo))
         θ = round(θ; digits=3)
#         @show θ
      end

      nparams[perm] .= nparams[perm] .+ β

      if STAGE==1
      #nvals,nvecs, = tsrcalc2(nparams,stg,cdo,ctrl.NFOLD,ctrl,nlist)
         @time nvals,nvecs = tsrcalc_1stg!(nvals,nvecs,nlist,σs,ctrl,nparams,stg,ℋ)
      elseif STAGE==2
         nvals,nvecs,tvcs, = twostg_calc2(nparams,stg,cdo,ctrl.NFOLD,ctrl,nlist)
      else
         println("Someday....")
      end
      nrms, nomc, = rmscalc(nvals,inds,ofreqs)
      nwrms = √(nomc' *W*nomc ./ length(nomc))
      #ρlm = lbmq_gain(β,λlm,jtw,H,omc,nomc)
      ρlm = lbmq_gain2(β, J,omc,nomc)
      check = (nrms-rms)/rms
#      @show check
      @show norm(β)
      println("ρlm = $(round(ρlm; digits=4)), nrms = $(round(nrms; digits=4)), wrms = $(round(nwrms; digits=4)), θ = $θ")
#      println("Δlm = $(round(Δlm; digits=4)), wrms = $(round(nwrms; digits=4))")
#      println("λlm = $(round(λlm; digits=4)), norm(β) = $(round(norm(β);digits=4))")
#      println()
      #println("norm(β ./ params) = $(norm(β ./ params[perm]))")
      #println("inv(H)*β = $(inv(H)*β)")
      #if (nrms<rms)&&(ρlm>1e-3)
      #stepcheck = (ρlm>1e-5)
      if (BOLD == 0)||wrms > 1e5
         stepcheck = ((nrms<rms)&&(ρlm>1e-4))
      else
      stepcheck = ((nrms<rms)&&(ρlm>1e-2)) || ((nrms*(1-θ)^BOLD)<0.2*lrms)
      end
      #stepcheck = check > 1e-8
      @show stepcheck
      #if ((nrms*(1-θ)^BOLD)< lrms)&&BOLD>0|| ((nrms<rms)&&(ρlm>1e-3))# || bad > 2
      if stepcheck
      #@show θ
      #@show round.(bf, sigdigits=5)
      #@show round.(bf./nparams[perm], sigdigits=5)
         counter += 1
         rms = nrms
         lrms = min(lrms,rms)
         omc = nomc
         tΔ = round.((nvals .- vals)./sum(βf), sigdigits=5)
         vals .= nvals
         params[perm] .+= β
         vecs .= nvecs
         wrms = √(omc' *W*omc ./ length(omc))
         #println(params)
         #@time J = build_jcbn!(J,cdo,inds,S,ctrl,vecs,params,perm,scales)
      if STAGE==1
      #@time J = build_jcbn2!(J,cdo,nlist,inds,ctrl,nvecs,nparams,perm,scales,stg)
         printstyled("updaing jacobian\n",color=:green)
         @time build_jcbn!(J,inds,nlist,ℋ,ctrl,perm,vecs,params,scales,stg)
      elseif STAGE==2
      @time J = build_jcbn_2stg!(J,cdo,nlist,inds,ctrl,nvecs,nparams,perm,scales,stg,tvcs)
      else
         println("Someday....")
      end      
         #J, w = linereject(J,W,omc,uncs,ctrl.REJECT)
         #build_hess!(H,jtw,J,W)
   build_hess!(H,jtw,J,W)
         βo = β
            #println(diag(H))
         paramarray[:,counter+1] = params
         #sρlm = (@sprintf("%0.4f", ρlm))
         srms = (@sprintf("%0.4f", rms))
         slλ = (@sprintf("%0.4f", log10(λlm)))
         scounter = lpad(counter,3)
         println("After $scounter iterations, RMS = $srms, log₁₀(λ) = $slλ")
         #iterationwriter(molnam,paramarray,rms,counter,λlm,sum(βf,dims=2),perm)
#         @show wrms
#         @show params
         μlm /= 10.0
         Δlm *= min(1.1,5.0)
         GEO=false#true
      else
         GEO = false
         #GEO = GEO != true #inverts GEO
         μlm = max(4.0*μlm,1.0E-24)
         Δlm = max(Δlm*.8,1.0E-24)
         #λlm = max(10.0*λlm,1.0E-24)
      end
   converged,endpoint = fincheck!(converged,endpoint,rms,βf,
         λlm,goal,abs(check),ϵ0,ϵ1,ϵ2,counter,LIMIT,params[perm],-2jtw*omc)
   end#while
   frms, fomc, fcfrqs = rmscalc(vals, inds, ofreqs)
   puncs = zeros(size(params))
   puncs[perm] = paramunc(H,W,perm,omc)
   covarmat = covarr2(H,omc)#covarr(correl(H),puncs)
   #@show H
   #@show correl(H)
   #params[1:15] .= paramrecov(params[1:15])
   #uncs[1:15] .= uncrecov(uncs[1:15],params[1:15])
   params[1:18], puncs[1:18] = fullrecov(params[1:18],puncs[1:18],ctrl.Irrep)
   slλ = (@sprintf("%0.4f", log10(λlm)))
   outputfinal(molnam,ctrl,frms,counter,slλ,puncs,params,endpoint)
   triangleprint(correl(H);io="$molnam.out")
   if occursin("J", ctrl.RUNmode)&&STAGE==2
      out_jcbn_2stg(molnam,ops,nlist,inds,ctrl,vecs,params,scales,stg,tvcs)
   end
   if ctrl.overwrite==true
      println("Writing new input file at $molnam.inp. Previous file has moved to $molnam","1.inp")
      inpwriter(molnam, params, scales)
   end
   return params, covarmat, fomc, fcfrqs, vals
end
