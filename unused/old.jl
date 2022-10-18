"""
This collects functions that were previously used in the code but I seem to have
replaced or otherwise abandonded. Functions that I have future plans for are
included with their related functions
"""

function Htsr0N(pr,j,s,n,mcalc,σ)
   #array = NFOLD .* collect(Float64,-mcalc:mcalc) .+ σ
   marray = msbuilder(Float64,mcalc,σ,NFOLD)
   karray = collect(Float64,-n:n)
   if n == zero(n)
      mat = spzeros(Float64,length(marray),length(marray))#diagm(0=>ondiags)
   else
      mat = pr[12]*θ(j,n,s)*0.0
      mat .*= karray
      mat = kron(marray,mat)
      mat = spdiagm(0=>mat)
   end
   return mat
end
function Htsr1N(pr,j::Float64,s::Float64,nl,mcalc,σ)
   karray = collect(Float64,-nl:nl)
   #array = NFOLD .* collect(Float64,-mcalc:mcalc) .+ σ
   marray = msbuilder(Float64,mcalc,σ,NFOLD)
   mat = spzeros(Float64,0,0)
   for i in 1:length(marray)
      p1 = ϕ(j,nl+1.0,s)*pr[12]*marray[i]*0.0
      p1 .*= sqrt.( (nl+1.0)^2 .+ karray .^2)
      part = spdiagm((2*Int(nl)+1),(2*Int(nl)+3), 1=>p1)
      mat = cat(mat,part,dims=(1,2))
   end
   return mat
end

function Htsr(pr,J,S,mcalc,σ)#OLD
   md = 2*mcalc + 1
   ns, nd, ni, jd = srprep(J,S,md)
#   ni = ni
   ni[1,1] = 1 #investigate why this is needed should be definitional to srprep but isn't???
   out = spzeros(Float64,md*jd,md*jd)
   out[1:ni[1,2],1:ni[1,2]] = kron(eye(md),Hrot(pr,ns[1]) + Hspi0N(pr,J,S,ns[1]))
   out[1:ni[1,2],1:ni[1,2]] += Htor(pr,mcalc,ns[1],σ) + Htsr0N(pr,J,S,ns[1],mcalc,σ)
   for i in 2:length(ns)
   n = ns[i]
   n1part = kron(eye(md),Hspi1N(pr,J,S,n-1.0)) + Htsr1N(pr,J,S,n-1.0,mcalc,σ)
   @inbounds out[ni[i-1,1]:ni[i-1,2],   ni[i,1]:ni[i,2]] = n1part
   @inbounds out[   ni[i,1]:ni[i,2],   ni[i,1]:ni[i,2]] = kron(eye(md),
               Hrot(pr,n) + Hspi0N(pr,J,S,n)) + Htor(pr,mcalc,n,σ) + Htsr0N(pr,J,S,n,mcalc,σ)
   @inbounds out[   ni[i,1]:ni[i,2],ni[i-1,1]:ni[i-1,2]] = transpose(n1part)
   end
   return out
end

function lbmq_step2(jcbn, weights, omc, λ, perm)
"""
This should be the Levenberg-Marquadt step. This solves (JᵗWJ+λI)Δβ = (JᵗW)Δy
   for Δβ. Where J is the Jacobian, W is the weights, λ is the
   Levenberg-Marquadt parameter, and Δy is the omcs. This returns the LVMQ step
   (βlm), ST step (βsd), and the parameter t.
"""
   jcbn = jcbn[:,perm]
   jtw = transpose(jcbn)*weights
   βlm = zeros(size(perm))
   jtj = jtw*jcbn
   #A = Hermitian(jtj + λ*Diagonal(jtj))
   A = Hermitian(jtj + λ*I)
   A = factorize(Symmetric(A))
   X = -jtw*omc
   βlm = ldiv!(βlm, A, X)
   βsd = -transpose(jcbn)*omc
   t = norm(βsd)^2/(norm(jcbn*βsd)^2)
   return βlm, βsd, t, X
end
function lbmq_opt(nlist,ofreqs,uncs,inds,params,scales,λ)
   vals,vecs = limeigcalc(nlist, inds, params)
   rms, omc = rmscalc(vals, inds, ofreqs)
   #println(omc)
   perm,n = findnz(sparse(scales))
   println("Initial RMS = $rms")
   counter = 0
   goal = sum(uncs)/length(uncs)
   newparams = copy(params)
   weights = diagm(0=>(uncs .^ -1))
   #weights = diagm(0=>ones(size(uncs)))
   converged = false
   THRESHOLD = 1.0E-8
   RHOTHRES = -1.0E-6
   LIMIT = 20
   λ0 = λ
   Δₖ = 1.0
   Δ0ₖ = Δₖ
   λ0 = λ
   while (converged==false)
      #println("Building Jacobians")
      jcbn = build_jcbn(inds,vecs,params)
      lgscls = 10 .^ (floor.(log10.(abs.(params[perm] ./maximum(params[perm])))))
      #lgscls = ones(size(lgscls))
      #println("Performing Levenberg-Marquadt Step")
      if true
         adjst,g = lbmq_step(jcbn,weights,omc,λ,perm) #.* lgscls
         normadjst = abs(norm(adjst))
         #if normadjst > Δₖ
         #   adjst = adjst .* (Δₖ/normadjst)
         #end
      else
         println("dogleg")
         βlm, βsd, t, g = lbmq_step2(jcbn,weights,omc,λ,perm)
         adjst = dogleg(βlm,βsd,t,Δₖ)
      end
      adjst .*= scales[perm]
      #back up parameters
      newparams[perm] = params[perm] .+ adjst
      #recalculate RMS
      vals,nvecs = limeigcalc(nlist, inds, newparams)
      nrms, nomc = rmscalc(vals, inds, ofreqs)
      ρlm = lbmq_gain(adjst,λ,g,rms,nrms)
      println(ρlm)
      #println(adjst[1],"   ", adjst[end])
      check = abs(nrms-rms)/rms
      counter += 1
      if ρlm > 0.0#-1.0E-6 #nrms ≤ rms
         #accept step and decrease λ
         params = newparams
         rms = nrms
         omc = nomc
         vecs = nvecs
         Δₖ *= 1.5
         λ = λ/3.0 #max(1/3,1-(2*ρlm-1)^3)
         #νlm = 2.0
      #elseif (ρlm > RHOTHRES)&&(ρlm < 0.0)
      #   Δₖ = Δ0ₖ
      #   λ = λ0
      #   counter -= 4
      else #nrms > rms
         #reject step due to RMS increase
         λ = λ*2.0
         #λ = min(λ,1.0E+12)
         Δₖ *= 0.9
         #Δₖ = max(Δₖ,0.00001)
         #params[perm] = params[perm] .+ adjst
         #rms = nrms
      end #ρlm if
      srms = (@sprintf("%0.4f", rms))
      slλ = (@sprintf("%0.4f", log10(λ)))
      sΔ = (@sprintf("%0.6f", Δₖ))
      scounter = lpad(counter,3)
      println("After $scounter interations, RMS = $srms, log₁₀(λ) = $slλ, Δₖ = $sΔ")
      #println(check)
      if (check < THRESHOLD)#||(rms ≤ goal)#&&(counter > 1)
         println("A miracle has come to pass. The fit has converged")
         break
      elseif counter ≥ LIMIT
         println("Alas, the iteration count has exceeded the limit")
         #println(omc)
         break
      else
         #write update to file
      end #check if
   end#converged while
   #println(omc)
   return params, vals
end


function cost(vals,inds,ofreqs,weights)
   cfreqs = zeros(size(ofreqs))
   for i in 1:size(cfreqs)[1]
      cfreqs[i] = vals[inds[i,3],inds[i,2]+1] - vals[inds[i,6],inds[i,5]+1]
   end
   omc = ofreqs - cfreqs
   cst = sum((omc ./weights) .^2)
   rms = sqrt(cst/length(omc))
   return cst
end
function build_jcbn(inds,vecs,params)
"""
This builds the Jacobian based on the Hellmann–Feynman theorem.
"""
   jcbn = zeros(Float64,size(inds)[1],length(params))
   Threads.@threads for a in 1:size(inds)[1]
      ju = 0.5*inds[a,1]
      jl = 0.5*inds[a,4]
      σu = inds[a,2]
      σl = inds[a,5]
      vecu = vecs[1:Int((2*S+1)*(2*ju+1)*(2*mcalc+1)),inds[a,3],σu+1]
      vecl = vecs[1:Int((2*S+1)*(2*jl+1)*(2*mcalc+1)),inds[a,6],σl+1]
      for b in 1:length(params)
         jcbn[a,b] = anaderiv(jl,S,σl,vecl,params,b) - anaderiv(ju,S,σu,vecu,params,b)
      end
   end
   return jcbn
end
function lbmq_step(jcbn, weights, omc, λ, perm)
"""
This should be the Levenberg-Marquadt step. This solves (JᵗWJ+λI)Δβ = (JᵗW)Δy
   for Δβ. Where J is the Jacobian, W is the weights, λ is the
   Levenberg-Marquadt parameter, and Δy is the omcs. This returns the step, Δβ.
"""
   jcbn = jcbn[:,perm]
   jtw = transpose(jcbn)*weights
   β = zeros(size(perm))
   jtj = jtw*jcbn
   A = Hermitian(jtj + λ*Diagonal(jtj))
   A = factorize(Symmetric(A))
   X = jtw*omc
   β = ldiv!(β, A, -X)
   return β,X
end

function qngen(n,m,σ)
   #depreciated
   nd = 2*n+1
   md = 2*m+1
   narray = fill(n,nd*md)
   karray = kron(ones(md),collect(-n:n))
   marray = kron(NFOLD .* collect(-m:m) .+ σ,ones(nd))
   σarray = fill(σ,nd*md)
   out = hcat(narray,karray,marray,σarray)
end
function qngen(j,s,m,σ)
   #depreciated
   nlist = Δlist(j,s)
   out = zeros(0,4)
   for i in 1:length(nlist)
      out = vcat(out,qngen(nlist[i],m,σ))
   end
   jlist = fill(j,size(out)[1])
   out = hcat(jlist,out)
   return out
end

function vecpadder(ns,degns,offst,nm,vtmi,vtc,vecs)
   #depreciated
   partial = Array{Float64}(undef,0,sum(degns))
   for i in 1:length(ns)
      pad = (nm - ns[i])
      zpad = zeros(Float64, pad, sum(degns))
      for v in 0:(vtc)
         temp = vecs[offst[i]+v*degns[i]+1:offst[i]+(v+1)*degns[i],:]
         temp = vcat(zpad,temp,zpad)
         partial = vcat(partial,temp)
      end
   end
   return partial
end

function greaterof(x,y)
"""
I don't remember why I needed this but it returns the greater of two input values.
"""
   if x > y
      return x
   else
      return y
   end
end

function intmat(jb,nb,jk,nk,s,k)
   sqn = length(k)
   mat = zeros(Float64,sqn,sqn)
   for x in 1:sqn
   for y in x:sqn
      @inbounds mat[x,y] = intelem(jb,nb[y],k[y],s,jk,nk[x],k[x])
   end
   end
   return Symmetric(mat)
end
function intbuild(jmax,mcalc,jb,jk,s)
   ns = Int(jmax + s)
   nbs = Δlist(jb,s)
   nks = Δlist(jk,s)
   nbarray = kron(ones((2*mcalc+1)*(2*ns+1)),nbs[1])
   nkarray = kron(ones((2*mcalc+1)*(2*ns+1)),nks[1])
   karray = kron(ones(2*mcalc+1),collect(-ns[1]:ns[1]))
   for i in 2:length(nbs)
      nbarray = vcat(nbarray,kron(ones((2*mcalc+1)*(2*ns+1)),nbs[i]))
      nkarray = vcat(nkarray,kron(ones((2*mcalc+1)*(2*ns+1)),nks[i]))
      karray = vcat(karray,kron(ones(2*mcalc+1),collect(-ns:ns)))
   end
   mat = intmat(jb,nbarray,jk,nkarray,s,karray)
   return mat
end
function intcalc(jmax,mcalc,s)
   sjmd = (2*s+1)*(2*jmax+1)*(2*mcalc+1)
   jind = convert(Int,jmax+s)
   μmat = zeros(Float64,sjmd,sjmd,jind,jind)
   for x in 1:jind
   for y in x:jind
      μmat[:,:,x,y] = intbuild(jmax,y-s,x-s,s)
      μmat[:,:,y,x] = μmat[:,:,x,y]
   end
   end
   return μmat
end



################################################################################
############                   unused at this time                  ############
################################################################################
function solvcalc(rprms,jlist,linds,ofreqs)
"""
This is the objective function for interfacing with other packages like Optim.jl
   It takes the parameters, the j σ pairs, line indices, and observed freqs.
   It then returns the rms. It can structured as x ->
   solvcalc(x,jlist,linds,ofreqs) for nice interfacing.
"""
   vals,vecs = limeigcalc(jlist, linds, rprms)
   rms, omc = rmscalc(vals, linds, ofreqs)
   println(rms)
   return rms
end
function solvcalc0(rprms,jlist,linds,ofreqs)
"""
This is the objective function for interfacing with other packages like Optim.jl
   It takes the parameters, the j σ pairs, line indices, and observed freqs.
   It then returns the rms. It can structured as x ->
   solvcalc(x,jlist,linds,ofreqs) for nice interfacing. It also prints the
   line indices, eigenvalues, and omcs for troubleshooting
"""
   vals,vecs = limeigcalc(jlist, linds, rprms)
   rms, omc = rmscalc(vals, linds, ofreqs)
   println(linds)
   println(vals)
   println(omc)
   return rms
end
function jcbn!(G,rps)
   G = build_jcbn(linds,vecs,rps)
   return G
end

function dHdA(j,s,σ,vec,rp)
   drp = zeros(Float64, length(rp)+1)
   A = rp[1]
   F = rp[5]
   #dAeff/dAp
   drp[1] = (F^2) / ((F - A)^2)
   #dFr/dAp
   drp[5] = (F^2) / ((F - A)^2)
   #dFρ/dAp
   drp[6] = (F^2) / ((F - A)^2)
   if σ==0
      U = ur(j,s,mcalc)
      U *= ut(mcalc,j,s)
   else
      mat = Matrix(Htsr(drp,j,s,mcalc,σ))
   end
   out = transpose(vec)*mat*vec
   return out[1]
end
function dHdF(j,s,σ,vec,rp)
   drp = zeros(Float64, length(rp)+1)
   A = rp[1]
   F = rp[5]
   #dAeff/dF
   drp[1] = -(A^2) / ((F - A)^2)
   #dFr/dF
   drp[5] = (F^2 - 2.0A*F) / ((F - A)^2)
   #dFρ/dF
   drp[6] = -(A^2) / ((F - A)^2)
   if σ==0
      U = ur(j,s,mcalc)
      U *= ut(mcalc,j,s)
   else
      mat = Matrix(Htsr(drp,j,s,mcalc,σ))
   end
   out = transpose(vec)*mat*vec
   return out[1]
end
function build_jcbn2(inds,vecs,params)
"""
This builds the Jacobian based on the Hellmann–Feynman theorem. This is a
   modified version that uses particular expressions for F and A due to their
   direct couplings through ρ
"""
   #fitprm = params[perm]
   jcbn = zeros(Float64,size(inds)[1],length(params))
   for a in 1:size(inds)[1]
      ju = 0.5*inds[a,1]
      jl = 0.5*inds[a,4]
      σu = inds[a,2]
      σl = inds[a,5]
      vecu = vecs[1:Int((2*S+1)*(2*ju+1)*(2*mcalc+1)),inds[a,3],σu+1]
      vecl = vecs[1:Int((2*S+1)*(2*jl+1)*(2*mcalc+1)),inds[a,6],σl+1]
      jcbn[1,b] = dHdA(ju,S,σu,vecu,params,b) - dHdA(jl,S,σl,vecl,params,b)
      jcbn[2,b] = anaderiv(ju,S,σu,vecu,params,b) - anaderiv(jl,S,σl,vecl,params,b)
      jcbn[3,b] = anaderiv(ju,S,σu,vecu,params,b) - anaderiv(jl,S,σl,vecl,params,b)
      jcbn[4,b] = anaderiv(ju,S,σu,vecu,params,b) - anaderiv(jl,S,σl,vecl,params,b)
      jcbn[5,b] = dHdF(ju,S,σu,vecu,params,b) - dHdF(jl,S,σl,vecl,params,b)
      for b in 6:length(params)
         jcbn[a,b] = anaderiv(ju,S,σu,vecu,params,b) - anaderiv(jl,S,σl,vecl,params,b)
      end
   end
   return jcbn
end

function westerfit_optim()
"""
This version of westerfit uses the Optim package currently with the Accelerated
   Gradient Descent. It is profoundly slow but does seem to approach convergence.
   Best to allow it to run bit by bit. It is programed to only run for about
   600 seconds at a time. It will then print the new parameters so it can be
   restarted using those values.
"""
   x0 = parameters
   #read the lines
#   lines = readdlm("$molnam.lne", ',', Float64)
   lines = testlines
   #determine the states
   linds, ofreqs, uncs = lineprep(lines)
   #println(linds)
   goal = sum(lines[:,end])/length(lines[:,end])
   jlist = jlister(linds)
   global jlist=jlist
   global linds = linds
   global ofreqs = ofreqs
   global mmax=0
   global mcalc=0
   global S=0.5
   global nmax= S + 0.5*maximum(jlist[:,1])
#   scales = vcat(ones(3),zeros(1))
#   println("Beginning optimization")
   init_rms = solvcalc0(x0,jlist,linds,ofreqs)
   println("Initial RMS = $init_rms")
   lower = [82755.0, -1.0E-37, -1.0E-37, -1.0E+6, 1.0E-37, -1.0E+09, -1.0E+09,
             254.3617999, -1.0E+09, -1.0E+09, -1.0E-37, -1.0E-37]
   upper = [82758.0, 1.0E+09, 1.0E+09, 1.0E+6, 1.0E+12, 1.0E+09, 1.0E+09, 254.3618001,
             1.0E+09, 1.0E+09, 1.0E-37, 1.0E-37]
#   res = optimize(solvcalc, tsrparams, ParticleSwarm(; n_particles=12))
#   inner_optimizer = GradientDescent(linesearch=LineSearches.BackTracking(order=3))
#   res = optimize(solvcalc, lower, upper, tsrparams, Fminbox(inner_optimizer))
#         Optim.Options(iterations=100))
   res = Optim.optimize(x->solvcalc(x,jlist,linds,ofreqs), x0, AcceleratedGradientDescent(),
   Optim.Options(time_limit=600))#,
   println(res)
   println("Final RMS = ", Optim.minimum(res), " MHz")
   rotparams = Optim.minimizer(res)
   println("New Parameter Vector:")
   println(rotparams)
#   println("New Energy levels")
#   for j in 1:nmax
#      qns, vals, vecs = tsrdiag(tsrparams,j,s,nmax,mcalc,0,σ)
#      println(vals)
#   end
   #write output file
end

function xtxsolve(A)
   x = zeros(size(A))
   a,b = eigen(A)
   x = √(diagm(a))*transpose(b)
   return x
end
