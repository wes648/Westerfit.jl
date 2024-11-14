function qfreqs!(cfreqs,vals,inds)
   @threads for i in 1:size(cfreqs,1)
      cfreqs[i] = vals[inds[i,3],inds[i,2]+1] - vals[inds[i,6],inds[i,5]+1]
   end
   return cfreqs
end

function outputstart(molnam,λ,rms)
   io = open("$molnam.out", "a")
   println(io,"Initial RMS = $rms MHz")
   println(io,"Initial λ = $λ")
   println(io,"")
   println(io,"-------------------------------------")
   println(io,"")
   close(io)
end

function iterprint(rms,λlm,counter)
   srms = (@sprintf("%0.4f", rms))
   slλ = (@sprintf("%0.4f", log10(λlm)))
   scounter = lpad(counter,3)
   println("After $scounter iterations, RMS = $srms, log₁₀(λ) = $slλ")
end

function lbmq_approx!(β,H,J,jtw,omc,λ,nparams,scls,perm)
   A = Hermitian(H + λ*Diagonal(H))
   while isposdef(A)==false #this could be tidier
      λ = max(2.0*λ,1.0E-24)
      #println(λ)
      A = Hermitian(H + λ*Diagonal(H))
   end
   #A = Hermitian(λ*Diagonal(H))
   if isinf(λ)
      @warn "LB-MQ Matrix not pos-def!"
      println("Make sure you aren't trying to optimize a parameter with value of 0.0.")
      println("Or Wes fucked up the Hessian builder again...")
      println("Either way, This code is about to crash")
   end
   A = cholesky!(A)
   β .= ldiv!(β, A, jtw*omc) .* scls[perm]
   @show β' * jtw*omc > zero(β' * jtw*omc)
   nparams[perm] .+= β
   return β,λ,nparams
end

function lbmq_approx(ctrl,nlist,ofreqs,uncs,inds,params,scales,cdo,stg,molnam)
   vals,vecs, = tsrcalc2(params,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
   cfreqs = zero(ofreqs)
   cfreqs = qfreqs!(cfreqs,vals,inds)
   LIMIT = ctrl["maxiter"]
   paramarray = zeros(Float64, length(params), LIMIT+1)
   paramarray[:,1]=params
   #oparams = paramarray[:,1]

   rms, omc, = rmscalc(vals, inds, ofreqs)
   perm = permdeterm(scales,stg)
   @show perm
   println("Initial RMS = $rms")
   goal = BLAS.nrm2(uncs)/√length(uncs)*ctrl["goal"]
   W = Diagonal(1.0 ./ uncs)
   #RHOTHRES = -1.0E-6
   ϵ0 = 0.1E-6 #RMS check
   ϵ1 = 0.1E-6 #step size check
   μlm = ctrl["λlm0"]#(rms + rms^2)#*0.0
   λlm = λgen(μlm, rms) 
   #oλlm = λlm
   println("Initial λ = $λlm")
   #Δlm = nrm2(params[perm])/length(perm)
   counter = 0

   outputstart(molnam,λlm,rms)
   nparams = copy(params)
   βf = zeros(length(perm)) #step
   J = zeros(Float64,size(inds,1),length(perm)) #Jacobian
   jtw = zero(omc) #gradient = -2jtw*omc
   J = build_jcbn3!(J,cdo,nlist,inds,ctrl,vecs,params,perm,scales,stg)
   #@show size(J)
   #J, w, omc = linereject(J,W,omc,uncs,ctrl["REJECT"])
   H, jtw = build_hess(jtw,J[:,perm],W)
   #@show H
   if true ∈ isnan.(H)
      println("FUCKING FUCKING FUCK. NaN in Hessian")
   end
   endpoint = "not yet"
   converged=false
while converged==false
   #@show sum(-2jtw*omc)
   #@show sum(βf)

   λlm = λgen(μlm, rms) 
   #βf,λ,nparams = lbmq_approx!(βf,H,J,jtw,omc,λlm,params,scales,perm)
   βf = J[:,perm] \ omc
   nparams[perm] += 0.5βf
   cfreqs = tsrapprox(J,nparams,stg)
   nrms, nomc = rmscalc(cfreqs,ofreqs)
   #nomc = ofreqs .- cfreqs
   ρlm = lbmq_gain(βf,λlm,-2jtw*omc,H,omc,nomc)
   #@show ρlm
   #check = abs(nrms-rms)/rms
   check = abs(nrms-rms)/rms

   #if (ρlm > 1.0e-7)#||(nrms<rms)#1.0e-7
   #if nrms < rms
      rms = nrms
      params .= nparams
      if (mod(counter,6)==0)&&(counter > 0)
         println("New vector time!")
         vals,vecs, = tsrcalc2(nparams,stg,cdo,ctrl["NFOLD"],ctrl,nlist)
         cfreqs = qfreqs!(cfreqs,vals,inds)
         @time J = build_jcbn3!(J,cdo,nlist,inds,ctrl,vecs,params,perm,scales,stg)
         H, jtw = build_hess(jtw,J[:,perm],W)
      else
         y = -2*jtw*(nomc - omc)
         H += (y*y')/(y' *βf) - (H*βf*βf' *H')/(βf' *H*βf)
#         cfreqs = tsrapprox(J,params)
      end
      omc = nomc
      counter += 1
      paramarray[:,counter+1] = params
      #@show norm(jtw*omc)
      iterprint(rms,λlm,counter)
      iterationwriter(molnam,paramarray,rms,counter,λlm,βf,perm)
      μlm /= 30.0
   #else
   #   printstyled("grumble grumble\n"; color=:red)
   #   μlm = max(4.0*μlm,1.0E-24)
   #end
   converged,endpoint = fincheck!(converged,endpoint,rms,βf,λlm,goal,check,ϵ0,ϵ1,counter,LIMIT,params[perm])
   if converged
      @show βf
      @show -2jtw*omc
      @show eigvals(H)
   end
end#while
   #@show H
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
