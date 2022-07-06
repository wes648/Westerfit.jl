



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
