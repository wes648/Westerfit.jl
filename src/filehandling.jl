#####INPUTS

function paraminp(molnam::String)
   NFOLD = 3
   S = 0.0
   TK = 25.0
   mcalc = 0
   mmax = 0
   Nmax = 0
   A = 0.0, 0.0
   B = 0.0, 0.0
   C = 0.0, 0.0
   Dab = 0.0, 0.0
   F = 0.0, 0.0
   ρ = 0.0, 0.0
   V3 = 0.0, 0.0
   ϵzz = 0.0, 0.0
   ϵxx = 0.0, 0.0
   ϵyy = 0.0, 0.0
   ϵxz = 0.0, 0.0
   η = 0.0, 0.0
   χzz = 0.0, 0.0
   χxmy = 0.0, 0.0
   χxz = 0.0, 0.0
   ΔN = 0.0, 0.0
   ΔNK = 0.0, 0.0
   ΔK = 0.0, 0.0
   δN = 0.0, 0.0
   δK = 0.0, 0.0
   Fm = 0.0, 0.0
   V6 = 0.0, 0.0
   V3m = 0.0, 0.0
   ρm = 0.0, 0.0
   ρ3 = 0.0, 0.0
   FN = 0.0, 0.0
   FK = 0.0, 0.0
   Fbc = 0.0, 0.0
   Fab = 0.0, 0.0
   V3N = 0.0, 0.0
   V3K = 0.0, 0.0
   V3ab = 0.0, 0.0
   V3bc = 0.0, 0.0
   ρN = 0.0, 0.0
   ρK = 0.0, 0.0
   ρab = 0.0, 0.0
   ρbN = 0.0, 0.0
   ΔsN = 0.0, 0.0
   ΔsNK = 0.0, 0.0
   ΔsKN = 0.0, 0.0
   ΔsK = 0.0, 0.0
   δsN = 0.0, 0.0
   δsK = 0.0, 0.0
   ΦJ = 0.0, 0.0
   ΦJK = 0.0, 0.0
   ΦKJ = 0.0, 0.0
   ΦK = 0.0, 0.0
   ϕJ = 0.0, 0.0
   ϕJK = 0.0, 0.0
   ϕK = 0.0, 0.0
   μa = 1.0, 0.008
   μb = 0.5, 0.006
   μc = 0.0, -0.034
   eval.(Meta.parse.(readlines(pwd()*"/"*molnam*".inp"))) #doesn't work due to scope behavior of eval()
   println(Nmax)
   if NFOLD==zero(NFOLD)
      if mcalc != zero(mcalc)
         print("NFOLD is zero; setting mcalc to 0")
         mcalc = 0
      end
      if mmax != zero(mmax)
         print("NFOLD is zero; setting mmax to 0")
         mmax = 0
      end
   end
   BJ = 0.5*(B[1]+C[1])
   BK = A[1] - BJ
   Bp = 0.25*(B[1] - C[1])
   ao = -(ϵzz[1] + ϵyy[1] + ϵxx[1])/3.0
   a = -(2.0*ϵzz[1] - ϵyy[1] - ϵxx[1])/6.0
   d = -ϵxz[1]*0.5
   b = (ϵxx[1] - ϵyy[1])*0.5
   χ2 = √(1.0/6.0)*χxmy[1]
   χ1 = -√(2.0/3.0)*χxz[1]
   params = [BK; BJ; Bp; Dab[1]; F[1]; ρ[1]*F[1]; V3[1]; ao[1]; a[1]; b[1]; d[1]
         η[1]; χzz[1]; χ2; χ1; ΔN[1]; ΔNK[1]; ΔK[1]; δN[1]; δK[1]; Fm[1]; V6[1]
         V3m[1]; ρm[1]; ρ3[1]; FN[1]; FK[1]; Fbc[1]; Fab[1]; V3N[1]; V3K[1]
         V3ab[1]; V3bc[1]; ρN[1]; ρK[1]; ρab[1]; ρbN[1]; ΔsN[1]; ΔsNK[1]
         ΔsKN[1]; ΔsK[1]; δsN[1]; δsK[1]; ΦJ[1]; ΦJK[1]; ΦKJ[1]; ΦK[1]; ϕJ[1]; ϕJK[1]; ϕK[1]]
   scales = [A[2]; B[2]; C[2]; Dab[2]; F[2]; ρ[2]*F[2]; V3[2]; ϵzz[2]; ϵxx[2]; ϵyy[2]; ϵxz[2];
           η[2]; χzz[2]; χxmy[2]; χxz[2]; ΔN[2]; ΔNK[2]; ΔK[2]; δN[2]; δK[2];
           Fm[2]; V6[2]; V3m[2]; ρm[2]; ρ3[2]; FN[2]; FK[2];
           Fbc[2]; Fab[2]; V3N[2]; V3K[2]; V3ab[2]; V3bc[2]; ρN[2]; ρK[2]; ρab[2]; ρbN[2]
           ΔsN[2]; ΔsNK[2]; ΔsKN[2]; ΔsK[2]; δsN[2]; δsK[2];
           ΦJ[2]; ΦJK[2]; ΦKJ[2]; ΦK[2]; ϕJ[2]; ϕJK[2]; ϕK[2]]
   μs = [μa[1] μa[2]; μb[1] μb[2]; 0.0 μc[2]]
   println(A)
   return params, scales, μs, Nmax, S, NFOLD, mcalc, mmax
end


function lineprep(lns,s,mcalc)
   #converts the input file into a more code friendly format
   #           1  2  3   4   5  6  7  8  9   10 11 12  13   14
   #input  = [ju nu kau kcu mu σu jl nl kal kcl ml σl freq unc]
   #           1  2   3   4  5   6
   #output = [ju σu indu jl σl indl]
   qunus = lns[:,1:12]
   freqs = lns[:,13]
   uncs = lns[:,end]
   inds = zeros(Int,size(lns)[1],6)
   inds[:,1] = Int.(2 .* qunus[:,1])
   inds[:,2] = Int.(qunus[:,6])
   inds[:,3] = qn2ind.(mcalc,qunus[:,5],qunus[:,1],s,qunus[:,2],qunus[:,3],qunus[:,4])
   inds[:,4] = Int.(2 .* qunus[:,7])
   inds[:,5] = Int.(qunus[:,12])
   inds[:,6] = qn2ind.(mcalc,qunus[:,11],qunus[:,7],s,qunus[:,8],qunus[:,9],qunus[:,10])
   #inds = vcat(inds[:,1:2], inds[:,3:4])
   return inds, freqs, uncs
end


#####OUTPUTS
function EngWriter(energies,qunus,mmax,sigma)
"""
Outputs energy levels with state assignments to a csv-like file
"""
   qunus[:,1] .*= 2
   qunus = convert.(Int,qunus)
   c = 29979.2458
   len = size(energies)[1]
   offsets = len .* collect(0:(2*mmax))
   out = fill("0",len)
   for i in 1:len
   #for m in 1:(2*mmax+1)
      energy = energies[i]/c
      #0.10f is what BELGI uses, 0.6f is for spcat
      part = lpad(@sprintf("%0.10f", energy), 16)
      part = string(part,",", lpad(@sprintf("%0.1f", qunus[i,1]/2), 7))
      out[i] = string(part,",", lpad(qunus[i,2],4),",",lpad(qunus[i,3],4),",",
       lpad(qunus[i,4],4),",", lpad(qunus[i,5],4),",", lpad(qunus[i,6],4))
#      lpad(qunus[i,5,1],4), lpad(sigma,4))
   #end
   end
   if sigma==0
      io = open("Astates_$molnam.eng", "w") do io
         for i in out
            println(io, i)
         end
      end
   elseif sigma==1
      io = open("Estates_$molnam.eng", "w") do io
         for i in out
            println(io, i)
         end
      end
   else
      println("Sorry, not ready for this σ value yet")
   end
end

function TraWriterSPCAT(molnam,freqs, qunus) #emulates the cat file structure of SPCAT
   c = 29979.2458
   p = sortperm(freqs[:,1])
   freqs = freqs[p,:]
   qunus = qunus[p,:]
   out = fill("0",size(freqs)[1])
   for i in 1:size(freqs)[1]
      #freq
#      part = lpad(@sprintf("%0.4f", freqs[i,1]), 13)
      part = @sprintf("%13.4f", freqs[i,1])
      #error
      part = string(part,@sprintf("%8.4f", 0.00))
      #-log(Intensity)
      modint = log(freqs[i,2])#*.1
      #modint = -abs(rand(1)[1])
      part = string(part, @sprintf("%8.4f", modint))
      #Degrees of Rotational Freedom
      part = string(part, lpad(3,2))
      #E_lower
      modEl = freqs[i,3]/c
      part = string(part, @sprintf("%10.4f", modEl))
      #Upper State degeneracy
      part = string(part, lpad(1,3))
      #Tag
      part = string(part, lpad(0,7))
      #QNFMT
      part = string(part, lpad(1415,4))
      #J N Ka Kc σ vt is the order in the array
      #N Ka Kc v J is the order for SPCAT
      #qunus for upper
      part = string(part, lpad(qunus[i,2],2),lpad(qunus[i,3],2),
      #lpad(qunus[i,4],2), lpad(qunus[i,6],2), lpad(qunus[i,5],2))
      lpad(qunus[i,4],2), lpad(qunus[i,5],2), lpad(qunus[i,1],2), lpad(qunus[i,6],2))
      #qunus for lower
      out[i] = string(part, lpad(qunus[i,8],2),lpad(qunus[i,9],2),
      lpad(qunus[i,10],2), lpad(qunus[i,11],2),  lpad(qunus[i,7],2), lpad(qunus[i,12],2))
      #lpad(qunus[i,10],2), lpad(qunus[i,12],2), lpad(qunus[i,11],2))
   end
   io = open("$molnam.cat", "w") do io
      for i in out
         println(io, i)
      end
   end
   println("Transitions written to $molnam.cat!")
end

function TraWriter(molnam,freqs, qunus)
   c = 29979.2458
   p = sortperm(freqs[:,1])
   freqs = freqs[p,:]
   qunus = qunus[p,:]
   out = fill("0",size(freqs)[1])
   counter = 0
   for i in 1:size(freqs)[1]
      #J N Ka Kc sigma vt is the order in the array
      #qunus for upper
      ju = @sprintf("%2.1f", 0.5*qunus[i,1])
      part = string(" ", lpad(ju,5),",", lpad(qunus[i,2],3),",", lpad(qunus[i,3],3),",",
      lpad(qunus[i,4],3),",", lpad(qunus[i,5],3),",", lpad(qunus[i,6],3))
      #qunus for lower
      jl = @sprintf("%2.1f", 0.5*qunus[i,7])
      part = string(part,",", lpad(jl,5),",", lpad(qunus[i,8],3),",", lpad(qunus[i,9],3),
      ",", lpad(qunus[i,10],3),",", lpad(qunus[i,11],3),",", lpad(qunus[i,12],3))
      #freq
#      part = lpad(@sprintf("%0.4f", freqs[i,1]), 13)
      freq = @sprintf("%13.4f", freqs[i,1])
      #error
      part = string(part,",", freq,",", @sprintf("%8.4f", 0.02))
      #-log(Intensity)
      #modint = log(freqs[i,2])#*.1
      modint = freqs[i,2]#*.1
      part = string(part, @sprintf("%8.4f", modint))
      if modint ≥ 0.0001
         counter += 1
         out[counter] = part
      end
   end
   out = out[1:counter]
   io = open("$molnam.cat", "w") do io
      for i in out
         println(io, i)
      end
   end
   println("Transitions written to $molnam.cat!")
end
