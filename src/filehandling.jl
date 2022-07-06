
function lineprep(lns)
   #converts the input file into a more code friendly format
   #           1  2  3   4   5  6  7  8  9   10 11 12  13   14
   #input  = [ju nu kau kcu mu σu jl nl kal kcl ml σl freq unc]
   #           1  2   3   4  5   6
   #output = [ju σu indu jl σl indl]
   qunus = lns[:,1:12]
   freqs = lns[:,13]
   uncs = lns[:,end]
   inds = zeros(Int64,size(lns)[1],6)
   inds[:,1] = Int64.(2 .* qunus[:,1])
   inds[:,2] = Int64.(qunus[:,6])
   inds[:,3] = qn2ind.(qunus[:,1],0.5,qunus[:,2],qunus[:,3],qunus[:,4])
   inds[:,4] = Int64.(2 .* qunus[:,7])
   inds[:,5] = Int64.(qunus[:,12])
   inds[:,6] = qn2ind.(qunus[:,7],0.5,qunus[:,8],qunus[:,9],qunus[:,10])
   #inds = vcat(inds[:,1:2], inds[:,3:4])
   return inds, freqs, uncs
end


#####OUTPUTS
function EngWriter(energies,qunus,mmax,sigma)
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
      #J N Ka Kc sigma vt is the order in the array
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

function TraWriter(molnam,freqs, qunus) #emulates the cat file structure of SPCAT
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
