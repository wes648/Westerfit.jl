#reduced barrier height is s = 4 V / F nf^2 (follows from the first term of the 0.5cosnx taylor series)

"""
This handles the output files for the westerfit package. Sophie made most of this!
"""

function paramrecov(prd::Array{Float64},ℋ::Vector{Op})::Array{Float64}
   out = zeros(18)

   #first loop over nfs to correct A, B, ρs, F, Vn_1s

   if prd[13] != 0.0
      out[14] = prd[14]/(-2.0*prd[13])                #ρ
      out[15] = prd[15]/(-1*prd[13])                  #ρx
   else
      out[14] = 0.0
      out[15] = 0.0
   end
   out[ 1] = prd[1] + prd[2] - prd[13]*out[14]^2        #A
   out[ 2] = prd[2] + 2.0*prd[3] - prd[13]*out[15]^2    #B
   out[ 3] = prd[2] - 2.0*prd[3]                        #C
   out[ 4] = prd[4] #- 2*out[14]*out[15]*csl             #Dab
   out[ 5] = (-prd[5] + √2*prd[7])/√3                   #ϵzz
   out[ 6] = -(prd[5]/√3 + prd[7]/√6) + prd[9]          #ϵxx
   out[ 7] = -(prd[5]/√3 + prd[7]/√6) - prd[9]          #ϵyy
   out[ 8] = -prd[8]                                    #ezx
   out[ 9] = prd[9]                                   #χzz
   out[10] = -√(1.5)*prd[12]                           #χxz
   out[11] = √6*prd[11]                                #χxx-χyy

   #for i ∈ 12:length(prd)
   #   if occursin("F",ℋ[i-11].nam)
   #      prd[i] /= csl
   #   elseif
   #   else
   #      break
   #   end 
   #end
   out[12] = prd[12]#/csl                        #F (MHz)
   out[14] = 2.0*prd[14] #/ csl                  #V3
   out[15] = prd[15]                             #η
   return out
end


function paramrecov_old(prd::Array{Float64})::Array{Float64}
   out = zeros(15)
   if prd[12] != 0.0
      out[13] = prd[13]/(-2.0*prd[12])       #ρ
   else
      out[13] = 0.0
   end
   out[1] = prd[1] + prd[2] - prd[12]*out[13]^2  #A
   out[2] = prd[2] + 2.0*prd[3]                  #B
   out[3] = prd[2] - 2.0*prd[3]                  #C
   out[4] = prd[4]                               #Dab
   out[5] = (-prd[5] + √2*prd[6])/√3             #ϵzz
   out[6] = -(prd[5]/√3 + prd[6]/√6) + prd[8]    #ϵxx
   out[7] = -prd[7]                              #ϵxz
   out[8] = -(prd[5]/√3 + prd[6]/√6) - prd[8]    #ϵyy
   out[9] = prd[9]                               #χzz
   out[10] = -√(1.5)*prd[10]                     #χxz
   out[11] = √6*prd[11]                          #χxx-χyy
   out[12] = prd[12]#/csl                        #F (MHz)
   out[14] = 2.0*prd[14] #/ csl                  #V3
   out[15] = prd[15]                             #η
   return out
end

function paramrecov_full(prd::Array{Float64})::Array{Float64}
   out = zeros(18)
   if prd[13] != 0.0
      out[14] = prd[14]/(-2.0*prd[13])                #ρ
      out[15] = prd[15]/(-1*prd[13])                  #ρx
   else
      out[14] = 0.0
      out[15] = 0.0
   end
   out[1] = prd[1] + prd[2] - prd[13]*out[14]^2        #A
   out[2] = prd[2] + 2.0*prd[3] - prd[13]*out[15]^2    #B
   out[3] = prd[2] - 2.0*prd[3]                        #C
   out[4] = prd[4] #- 2*out[14]*out[15]*csl             #Dab
   out[5] = (-prd[5] + √2*prd[7])/√3                   #ϵzz
   out[6] = -(prd[5]/√3 + prd[7]/√6) + prd[9]          #ϵxx
   out[7] = -(prd[5]/√3 + prd[7]/√6) - prd[9]          #ϵyy
   out[8] = (prd[6]-prd[8])                            #ezx
   out[9] = -(prd[6]+prd[8])                           #ϵxz
   out[10] = prd[10]                                   #χzz
   out[11] = -√(1.5)*prd[11]                           #χxz
   out[12] = √6*prd[12]                                #χxx-χyy
   out[13] = prd[13]#/csl                               #F (MHz)
   out[16] = 2.0*prd[16] #/ csl                         #V3
   out[17] = prd[17]                                   #η
   out[18] = 2*prd[18]                                 #ηx
   return out
end

function uncrecov(unc,prd::Array{Float64})::Array{Float64}
   out = zeros(15)

   if prd[12] != 0.0
      out[13] = (0.5*unc[13]/prd[12])^2 +
                (0.5*prd[13]*unc[12]/prd[12]^2)^2      #σρ
   else
      out[13] = 0.0
   end
#   out[1] = unc[1]^2 + unc[2]^2 + (0.25unc[12]*prd[13]^2/prd[12]^2)^2 +
#            (0.5*prd[13]*unc[13]/prd[12])^2          #σA
   out[1] = unc[1]^2 + unc[2]^2 + (unc[12]*prd[13]^2)^2 +
            (2*prd[12]*prd[13]*out[13])^2               #σA
   out[2] = unc[2]^2 + 4.0*unc[3]^2                     #σB
   out[3] = unc[2]^2 + 4.0*unc[3]^2                     #σC
   out[4] = unc[4]^2                                    #σDab
   out[5] = (unc[5]^2 + 2*unc[6]^2)/3                   #σϵzz
   out[6] = (2*unc[5]^2 + unc[6]^2)/6 + unc[8]^2        #σϵxx
   out[7] = unc[7]^2                                    #σϵxz
   out[8] = (2*unc[5]^2 + unc[6]^2)/6 + unc[8]^2        #σϵyy
   out[9] = unc[9]^2                                    #σχzz
   out[10] = 1.5*unc[10]^2                              #σχxz
   out[11] = 6.0*unc[11]^2                              #σχxx-χyy
   out[12] = unc[12]^2 #/csl                            #σF
   out[14] = 4.0*unc[14]^2 #/ csl                       #σV3
   out[15] = unc[15]^2                                  #ση
   return sqrt.(out)
end


function uncrecov_full(unc,prd::Array{Float64})::Array{Float64}
   out = zeros(18)

   if prd[13] != 0.0
      out[14] = (0.5*unc[14]/prd[13])^2 +
                (0.5*prd[14]*unc[13]/prd[13]^2)^2      #σρz
      out[15] = (0.5*unc[15]/prd[13])^2 +
                (0.5*prd[15]*unc[13]/prd[13]^2)^2      #σρx
   else
      out[14] = 0.0
      out[15] = 0.0
   end
#   out[1] = unc[1]^2 + unc[2]^2 + (0.25unc[12]*prd[13]^2/prd[12]^2)^2 +
#            (0.5*prd[13]*unc[13]/prd[12])^2          #σA
   out[1] = unc[1]^2 + unc[2]^2 + (unc[13]*prd[14]^2)^2 +
            (2*prd[13]*prd[14]*out[14])^2                       #σA
   out[2] = unc[2]^2 + 4.0*unc[3]^2+ (unc[13]*prd[15]^2)^2 +
            (2*prd[13]*prd[15]*out[15])^2                       #σB
   out[3] = unc[2]^2 + 4.0*unc[3]^2                             #σC
   out[4] = unc[4]^2                                    #σDab
   out[5] = (unc[5]^2 + 2*unc[6]^2)/3                   #σϵzz
   out[6] = (2*unc[5]^2 + unc[6]^2)/6 + unc[9]^2        #σϵxx
   out[7] = (2*unc[5]^2 + unc[6]^2)/6 + unc[9]^2        #σϵyy
   out[8] = unc[6]^2 + unc[8]^2                         #σϵzx
   out[9] = unc[6]^2 + unc[8]^2                         #σϵxz
   out[10] = unc[10]^2                                  #σχzz
   out[11] = 1.5*unc[11]^2                              #σχxz
   out[12] = 6.0*unc[12]^2                              #σχxx-χyy
   out[13] = unc[13]^2 #/csl                            #σF
   out[16] = 4.0*unc[16]^2 #/ csl                       #σV3
   out[17] = unc[17]^2                                  #ση
   out[18] = 4.0*unc[18]^2                              #σηx
   return sqrt.(out)
end

function fullrecov(prd,unc,irrep)
   oprd = paramrecov_full(prd)
   ounc = uncrecov_full(unc,oprd)
   #oprd[1:3] = oprd[irrepswap(irrep)]
   #ounc[1:3] = ounc[irrepswap(irrep)]
   printstyled("wes removed irrepswap here\n",color=:light_cyan)
   oprd[[13,16]] ./= csl
   ounc[[13,16]] ./= csl
   return oprd, ounc
end

#####OUTPUTS
function reslin(line,omc,cfrq)
   part  = lpad(line[ 1],4)*";"
   part *= lpad(Int(line[ 2]),3)*";"
   part *= lpad(Int(line[ 3]),3)*";"
   part *= lpad(Int(line[ 4]),3)*";"
   part *= lpad(Int(line[ 5]),3)*";"
   part *= lpad(line[ 6],4)*";"
   part *= lpad(Int(line[ 7]),3)*";"
   part *= lpad(Int(line[ 8]),3)*";"
   part *= lpad(Int(line[ 9]),3)*";"
   part *= lpad(Int(line[10]),3)*";"
   part *= " "*lpad(@sprintf("%0.5f", line[11]), 16)*";"
   part *= " "*lpad(@sprintf("%0.6f", omc), 16)*";"
   part *= " "*lpad(@sprintf("%0.5f", cfrq), 16)
end

function reswritter(molnam,lines,omcs,cfrqs)
   len = size(omcs,1)
   out = fill("0", len)
   @simd for i in 1:len
      out[i] = reslin(lines[i,:],omcs[i],cfrqs[i])
   end
   io = open("$molnam.res", "w") do io
      for i in out
         println(io, i)
      end
   end
end

function englin(s,eng,qunl)
   if s==zero(s)
      part = lpad(qunl[2],4)*","
   else
     # part = lpad(qunl[1],4)*"/2,"
      part = lpad(qunl[1],4)*","
      part *= lpad(qunl[2],4)*","
   end
   part *= lpad(qunl[3],4)*","
   part *= lpad(qunl[4],4)*","
   part *= lpad(qunl[5],4)*","
   part *= lpad(qunl[6],4)*","
   part *= " "*lpad(@sprintf("%0.10f", eng), 16)
   return part
end
function englin(s,eng,qunl,pasz)
   if s==zero(s)
      part = lpad(qunl[2],4)*","
   else
     # part = lpad(qunl[1],4)*"/2,"
      part = lpad(qunl[1],4)*","
      part *= lpad(qunl[2],4)*","
   end
   part *= lpad(qunl[3],4)*","
   part *= lpad(qunl[4],4)*","
   part *= lpad(qunl[5],4)*","
   part *= lpad(qunl[6],4)*","
   part *= " "*lpad(@sprintf("%0.10f", eng), 16)
   part *= ", "*lpad(@sprintf("%0.10f", pasz), 16)
   return part
end
function engwriter(molnam,ctrl,energies,qunus)
"""
Outputs energy levels with state assignments to a csv-like file
"""
   c = 29979.2458
#   @show size(energies)
   eng = energies[:,1]
   qns = qunus[:,:,1]
   for sc in 2:σcount(ctrl.NFOLD[1])
      eng = vcat(eng,energies[:,sc])
      qns = vcat(qns,qunus[:,:,sc])
   end
   len = size(eng,1)
   out = fill("0",len)
   for i in 1:len
      energy = eng[i]/c
      #0.10f is what BELGI uses, 0.6f is for spcat
      out[i] = englin(ctrl.S,energy,qns[i,:])
   end
   io = open("$molnam.eng", "w") do io
      for i in out
         println(io, i)
      end
   end
end
function engwriter(molnam,ctrl,energies,qunus,pasz)
"""
Outputs energy levels with state assignments to a csv-like file
"""
   c = 29979.2458
   eng = energies[:,1]
   qns = qunus[:,:,1]
   psz = pasz[:,1]
   for sc in 2:σcount(ctrl.NFOLD)
      eng = vcat(eng,energies[:,sc])
      qns = vcat(qns,qunus[:,:,sc])
      psz = vcat(psz,pasz[:,sc])
   end
   len = size(eng,1)
   out = fill("0",len)
   for i in 1:len
      energy = eng[i]/c
      #0.10f is what BELGI uses, 0.6f is for spcat
      out[i] = englin(ctrl.S,energy,qns[i,:],pasz[i])
   end
   io = open("$molnam.eng", "w") do io
      for i in out
         println(io, i)
      end
   end
end

function linestrng(s,frql,qunl)
   if s==zero(s)
      part  = lpad(qunl[2],3)*","  #N
      part *= lpad(qunl[3],3)*","  #Ka
      part *= lpad(qunl[4],3)*","  #Kc
      part *= lpad(qunl[5],3)*","  #m
      part *= lpad(qunl[8],3)*","  #N
      part *= lpad(qunl[9],3)*","  #Ka
      part *= lpad(qunl[10],3)*"," #Kc
      part *= lpad(qunl[11],3)*"," #m
      part *= " "*@sprintf("%13.4f", frql[1])*","
      part *= @sprintf("%10.4f", frql[4])*","
      part *= @sprintf("%12.6f", frql[2])*","
      part *= @sprintf("%10.4f", frql[3])
   else
      part  = lpad(qunl[1]*0.5,4)*"," #J
      part *= lpad(qunl[2],3)*","     #N
      part *= lpad(qunl[3],3)*","     #Ka
      part *= lpad(qunl[4],3)*","     #Kc
      part *= lpad(qunl[5],3)*","     #m 
      part *= lpad(qunl[7]*0.5,4)*"," #J
      part *= lpad(qunl[8],3)*","     #N
      part *= lpad(qunl[9],3)*","     #Ka
      part *= lpad(qunl[10],3)*","    #Kc
      part *= lpad(qunl[11],3)*","    #m
      part *= " "*@sprintf("%13.4f", frql[1])*","
      part *= @sprintf("%10.4f", frql[4])*","
      part *= @sprintf("%12.6f", frql[2])*","
      part *= @sprintf("%10.4f", frql[3])
   end
   return part
end
function qnlinSPCAT(qunus,s)
      #J N Ka Kc m σ is the order in the array
      #N Ka Kc v J is the order for SPCAT
   if s==zero(s)
      part = lpad(1404,4)
      part *= @sprintf("%2i", qunus[2])           #N
      part *= @sprintf("%2i", qunus[3])           #Ka
      part *= @sprintf("%2i", qunus[4])           #Kc
      part *= @sprintf("%2i", (qunus[5]))      #m
      part *= lpad(qunus[8],6)           #N
      part *= @sprintf("%2i", qunus[9])           #Ka
      part *= @sprintf("%2i", qunus[10])          #Kc
      part *= @sprintf("%2i", (qunus[11]))     #m
   else
      part = lpad(1415,4)
      part *= @sprintf("%2i", qunus[2])           #N
      part *= @sprintf("%2i", qunus[3])           #Ka
      part *= @sprintf("%2i", qunus[4])           #Kc
      part *= @sprintf("%2i", abs(qunus[5]))      #m
      part *= @sprintf("%2i", ceil(Int,qunus[1])) #J
      part *= "     "
      part *= @sprintf("%2i", qunus[8])           #N
      part *= @sprintf("%2i", qunus[9])           #Ka
      part *= @sprintf("%2i", qunus[10])          #Kc
      part *= @sprintf("%2i", abs(qunus[11]))     #m
      part *= @sprintf("%2i", ceil(Int,qunus[7])) #J
   end
   return part
end
function TraWriterSPCAT(molnam,s,freqs, qunus) #emulates the cat file structure of SPCAT
   c = 29979.2458
   p = sortperm(freqs[:,1])
   freqs = freqs[p,:]
   qunus = qunus[p,:]
   out = fill("0",size(freqs,1))
   for i in 1:size(freqs,1)
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
      part = string(part, @sprintf("%2i",3))
      #E_lower
      modEl = freqs[i,3]#/c
      part = string(part, @sprintf("%10.4f", modEl))
      #Upper State degeneracy
      degen = Int(2*qunus[i,1]+1)
      part = string(part, @sprintf("%3i",degen))
      #Tag
      part = string(part, @sprintf("%7i",7))
      #=#QNFMT
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
      #lpad(qunus[i,10],2), lpad(qunus[i,12],2), lpad(qunus[i,11],2))=#
      part *= qnlinSPCAT(qunus[i,:],s)
      out[i] = part
   end
   io = open("$molnam.cat", "w") do io
      for i in out
         println(io, i)
      end
   end
   println("Transitions written to $molnam.cat!")
end
function TraWriter(molnam,s, freqs, qunus)
   p = sortperm(freqs[:,1])
   freqs = freqs[p,:]
   qunus = qunus[p,:]
   out = fill("0",size(freqs,1))
   for i in 1:size(freqs,1)
      out[i] = linestrng(s, freqs[i,:], qunus[i,:])
   end
   io = open("$molnam.sim", "w") do io
      for i in out
         println(io, i)
      end
   end
   println("Transitions written to $molnam.cat!")
end


"""START CODE FOR UNCERTAINTY PRINTER THINGY"""

function num_to_string(x,fmt="%.1f")
   Printf.format(Printf.Format(fmt), x)
end

function uncrformatter(values,unc)
   uncertainty_digits = 3
   uncr = round.(unc, sigdigits = uncertainty_digits)
   uncstr = zeros(Float64, length(uncr))
   uncstr = string.(uncstr)

   for i in 1:length(uncr)
      number = -1*floor(Int, log10(unc[i])) + uncertainty_digits - 1
      words = string("%0.",number,"f")
      uncstr[i] = num_to_string(uncr[i],words)
   end

   valstr = zeros(Float64, length(values))
   valstr = string.(valstr)

   for i in 1:length(values)
      number = -1*floor(Int, log10(unc[i])) + uncertainty_digits - 1
      words = string("%0.",number,"f")
      valstr[i] = num_to_string(values[i],words)
   end

   uncstr1 = Base.lstrip.(uncstr1, '0')
   uncstr1 = Base.lstrip.(uncstr1, '.')
   uncstr1 = Base.lstrip.(uncstr1, '0')

   valunc = string.(valstr, "(", uncstr1, ")")
   return valunc
end

function uncrformattersci(values,unc)
   uncertainty_digits = 2

   uncr = round.(unc, sigdigits = uncertainty_digits)
   uncstr = fill("0", length(uncr))

   for i in 1:length(uncr)
      if uncr[i] == 0.0
         uncstr[i] = string("Fixed")
      elseif uncr[i] ≥ 1.0 && typeof(uncr[i]) == Float64
         uncstr[i] = num_to_string(uncr[i],"%0.2e")
      else
         number = abs(floor(Int, log10(unc[i]))) + uncertainty_digits - 1
         words = string("%0.",number,"f")
         uncstr[i] = num_to_string(uncr[i],words)
      end
   end
   
   uncstr1 = Base.lstrip.(uncstr, '0')
   uncstr1 = Base.lstrip.(uncstr1, '.')
   uncstr1 = Base.lstrip.(uncstr1, '0')
   uncunstr = tryparse.(Float64, uncstr1)

   for i in 1:length(uncunstr)
     # if typeof(uncunstr[i]) == Float64 && 10.0 <= uncunstr[i] < 100.0
     #    uncunstr[i] *= 10
      if typeof(uncunstr[i]) == Float64 && 1.0 <= uncunstr[i] < 10.0
         uncunstr[i] *= 100
      end
   end
   uncstr1 = string.(uncunstr)
   uncstr1 = Base.rstrip.(uncstr1, '0')
   uncstr1 = Base.rstrip.(uncstr1, '.')
   for i in 1:length(uncstr1)
      if uncstr1[i] != "nothing" && parse(Float64, uncstr1[i]) > 10.0 
         uncstr1[i] = Base.rstrip(uncstr1[i], '0')
      end
   end

   valstr = fill("0", length(values))
   
   for i in 1:length(values)
      if uncr[i] == 0.0
         number = length(values[i]) + 4
      else
         number = abs(floor(Int, log10(uncr[i]))) + uncertainty_digits - 1
      end
      words = string("%0.",number,"f")
      valstr[i] = num_to_string(values[i], words)
   end

   for i in 1:length(values)
      if uncr[i] == 0.0
         number = length(values[i]) + 4
      else
         number = -1*floor(Int, log10(uncr[i])) + uncertainty_digits - 1
      end
      temp = round(values[i], digits = number)
      valstr[i] = num_to_string(temp, "%0.10e")
   end
   
   valhalf = chop.(valstr, head = 0, tail = 4)
   valhalf = Base.rstrip.(valhalf, '0')
   
   ehalf = fill("0", length(values))
   for i in 1:length(valstr)
         ehalf[i] = chop(valstr[i], head = (length(valstr[i]) -4), tail = 0)
   end
   
   valunc = fill("0", length(values))
   
   finalval = fill("0", length(values))
   for i in 1:length(valstr)
      finalval[i] = chop(valstr[i], head = (length(valstr[i]) - 3), tail = 0)
   end
   finalval = tryparse.(Float64, finalval)
   evalue = fill("0", length(values))
   
   for i in 1:length(values)
      if finalval[i] % 3 == 0
         evalue[i] = num_to_string(finalval[i], "%.0f")
         valunc[i] = string(valhalf[i], "(", uncstr1[i], ")e", evalue[i])
      elseif (finalval[i] - 1) % 3 == 0
         valhalfreal = parse(Float64, valhalf[i])
         valhalfreal = valhalfreal*10
         valhalfreal = round(valhalfreal, sigdigits = length(valhalf[i]))
         finalval[i] = finalval[i] - 1
         evalue[i] = num_to_string(finalval[i], "%.0f")
         valhalf[i] = string(valhalfreal)
         valunc[i] = string(valhalf[i], "(", uncstr1[i], ")e", evalue[i])
      else (finalval[i] + 1) % 3 == 0
         valhalfreal = parse(Float64, valhalf[i])
         valhalfreal = valhalfreal * 100
         valhalfreal = round(valhalfreal, sigdigits = length(valhalf[i]))
         finalval[i] = finalval[i] - 2
         evalue[i] = num_to_string(finalval[i], "%.0f")
         valhalf[i] = string(valhalfreal)
         valunc[i] = string(valhalf[i], "(", uncstr1[i], ")e", evalue[i])
      end
   end

   for i in 1:length(values)
      if valhalf[i] == "0."
         valhalf[i] *= "0"
      end
      if valhalf[i] == "-0."
         valhalf[i] *= "0"
      end
      if uncr[i] == 0.0
         valunc[i] = string(valhalf[i], "e", evalue[i], " (Fixed)")
      elseif abs(unc[i]) >= abs(values[i])
         #valshort = num_to_string(valhalf[i], "%0.3f")
         valshort = valhalf[i][1:3]
         valunc[i] = string(valshort,"e",evalue[i],"(Undetermined)")
      else
      end
   end
   return valunc
end
"""END CODE FOR UNCERTAINTY PRINTER THINGY"""

"""INSERT CODE FOR INPUT FILE WRITER"""

function iterativecopier(molnam::String)
   oldarray = ["8", "7", "6", "5", "4", "3", "2", "1", ""]
   newarray = ["9", "8", "7", "6", "5", "4", "3", "2", "1"]
   for i in 1:length(oldarray)
      filename = string(molnam, oldarray[i], ".inp")
      newfilename = string(molnam, newarray[i], ".inp")
      if isfile(filename) == true
         cp(filename, newfilename, force=true)
      end
   end
end

function inpwriter2(molnam::String,ctrl,values)
   iterativecopier(molnam)

   #read 2ndorder block as ;-dlm file of strings
   #replace block[:,2] with strings.(values[2ndorder])
   block = readdlm("$molnam.inp",';',String,
                   skipstart=ctrl.sobk[1])[1:14+4*length(ctrl.NFOLD),1:3]
   block[:,2] .= lpad(string.(values[1:14+4*length(ctrl.NFOLD)]), 30)
   writedlm("$molnam.inp",block,';',skipstart=ctrl.sobk[1])

   #read params block as ;-dlm file of strings
   #replace block[:,2] with strings.(values[higherorder])
   block = readdlm("$molnam.inp",';',String, skipstart=ctrl.opbk[1])
   block[:,2] .= lpad(string.(values[14+4*length(ctrl.NFOLD):end]), 30)
   writedlm("$molnam.inp",block,';',skipstart=ctrl.opbk[1])

end

"""END CODE FOR INPUT FILE WRITER"""

function secnam_gen(nfold::Int)::Vector{String}
   l = length(nfold)
   out =  [" BK", " BN", " B⨦", " Dab", " T⁰₀(ϵ)", " T¹₁(ϵ)"," T²₀(ϵ)",
          " T²₁(ϵ)"," T²₂(ϵ)", " T²₀(χ)"," T²₁(χ)"," T²₂(χ)"]
   if isone(l)
      out = vcat(out,[" F", " -2ρzF", " -ρxF", " V$(nfold[1])/2", " ηz", " ηx"])
   elseif l>1
      for i ∈ 1:l
         out = vcat(out,[" F_$i", " -2ρzF_$i", " -ρxF$i", " V$(nfold[i])_$i/2", " ηz_$i", " ηx_$i"])
      end
   else
      #out is out
   end
   return out
end

function outputinit2(molnam,params,scls,linecount,ctrl,ℋ,rms,λ)
   secnam = secnam_gen(ctrl.NFOLD)
   io = open("$molnam.out","w")
      println(io,molnam,"   @   ",time)
      println(io,"")
      println(io,"Control Parameters")
      for fname in fieldnames(typeof(ctrl))
         println(io, "$fname = $(getfield(ctrl,fname))")
      end
      println(io,"")
      println(io,"Initial Parameters")
      l = length(secnam)
      for i in 1:l
         println(io, secnam[i], ";", lpad(param[i],30),";", lpad(scls[i],8))
      end
      println(io,"")
      for i in 1:length(ℋ)
         if scls[i] > 0.0
            println(io," ",ℋ[i].nam,";", lpad(param[i+l],30),";", lpad(scls[i+l],8))
         end
      end
      println(io,"\n\nNumber of lines = $linecount\n")
      println(io,"Initial RMS = $rms MHz\n","Initial λ = $λ\n\n",repeat("-",37),"\n")
   close(io)
end
function outputstart(molnam,λ,rms)
   outstr = string("Initial RMS = $rms MHz\n","Initial λ = $λ\n\n",repeat("-",37),"\n\n")
   io = open("$molnam.out", "a")
      print(io, outsrt)
   close(io)
end

function iterationwriter2(molnam,params,rms,counter,λlm,βf,perm)
   outstr = string("\nAfter $scounter Iterations:\n\n", 
      "RMS = $(@sprintf("%0.4f", rms)) MHz, log₁₀(λ) = $(@sprintf("%0.4f", log10(λlm)))\n\n",
      lpad("Parameter",33),lpad("Change",11),"\n")

   changes = fill("Fixed",length(params))
   changes[perm] .= @sprintf.("%0.4f",βf)
   nams = secnam_gen(ctrl.nfold)
   for i ∈ length(changes)
      ln = 30 - length(nams[i])
      outstr *= string(" ", nams[i],"; ", lpad(param[i],ln), "; ",lpad(changes[i],10),"\n")
   end

   ind = argmax(βf)
   maxnam = nams[perm][ind]
   percent = @sprintf("%0.4f", 100*βf[ind]/params[perm][ind] )

   outstr *= string("\nMax Change; $maxnam, $percent%\n", repeat("-",37), "\n")
   
   io = open("$molnam.out", "a")
      print(io, outstr)
   close(io)
end

function outputfinal(molnam,ctrl,ℋ,stgsfrms,counter,slλ,puncs,params,endpoint)

   strlnctrl,strln2nd,strlnhigh,file = findstrinput(molnam)

   srms = (@sprintf("%0.4f", frms))
   scounter = lpad(counter,3)

   #secnam = ["A","B","C","Dab","ϵzz","ϵxx","ϵzx","ϵxz","ϵyy","χzz","χxz",
   #          "χxx-χyy","F","ρz", "ρx", "V$(ctrl.NFOLD)","ηz","ηx"]
   #highnamall = file[strlnhigh:end,1]
   #highstg= file[strlnhigh:end,end]
   #highnam = highnamall[highstg .!= 0.0]
   #fullnam = vcat(secnam, highnam)

   fullnam = vcat(secnam_gen(nfold), ℋ[:].nam)

   if counter==0
      #This is a super sloppy bug fix for a bug when zero iterations occur
      outstr = repeat("-",111)
   else
      outstr = ""
   end
   if endpoint == "converge"
      outstr *= "\n A miracle has come to pass. The fit has converged.\n"
   elseif endpoint == "RMS"
      outstr *= "\n The RMS has stopped decreasing. Hopefully it is low.\n"
   elseif endpoint == "step size"
      outstr *= "\n It would appear step size has converged.\n"
   elseif endpoint == "LMthresh"
      outstr *= "\n λlm exceeded threshold.\n"
      outstr *= " If you were using the turducken, try again without it.\n"
   elseif endpoint == "iter"
      outstr *= "\n Alas, the iteration count has exceeded the limit.\n"
   elseif endpoint == "grad"
      outstr *= "\n The gradient has converged! Hurray!\n"
   else
      outstr *= "\n The output writer seems to be having issues.\n"
   end

outstr *= string("\nAfter $(lpad(counter,3)) iterations:\n\n","Final RMS = $srms MHz\n",
   "Final log₁₀(λ) = $slλ\n\n",lpad("Parameter",34),lpad("Uncertainty",31),lpad("‰ Uncertainty",17),"\n")

   for i in 1:length(params)
      if iszero(puncs[i])#puncs[i] == 0.0
         pdiffi = "-"
      else
         pdiffi = abs(round(puncs[i]/params[i]*1000,digits=4))
      end
      ln = 30 - length(fullnam[i])
      ln2 = 30 - length(params[i])
      ln3 = 17 - length(pdiffi)
      outstr *= string(" ",fullnam[i], "; ", lpad(params[i], ln), "; ", lpad(puncs[i], ln2),
          "; ", lpad(pdiffi, ln3),"\n")
   end
   outstr *= "\n\n"

   io = open("$molnam.out", "a")

      try
         finalunc = uncrformattersci(params,puncs)
         formatted = fill("0",2,length(fullnam))
         println(io,"Journal Formatted Values")
         for i in 1:length(fullnam)
            ln = 30 - length(fullnam[i])
            println(io, string(fullnam[i], "; ", lpad(finalunc[i],ln)) )
         end
      catch
         println(io," Yikes! The uncertainty formatter failed")
         println(io," There are likely some ill-fit parameters")
         println(io," You'll have to manually typeset the uncertainties\n")
      end

   println(io,"\n")
   if ctrl.apology # ==true
      println(io," Again, sorry about the name...")
   end
   println(io,"\n")
   close(io)

   println("Output written to $molnam.out!")
end


function triangleprint(mat,nams;io=stdout,d=4,col=5)
   io ≠ stdout ? io = open(io, "a") : io=stdout
   println(io,"  Correlation matrix:")
   l = size(mat,1)
   blocks = ceil(Int,l/col)
   for j in 1:blocks
   start = col*(j-1)+1
   stop = min(l,col*j)
   println(io,"\n"*prod(fill(" ",d))prod(lpad.(nams[start:stop],2d+2))*"\n")
   for i in start:l
      stop = min(i,col*j)
      part = lpad(nams[i],4)*prod(lpad.(round.(mat[start:stop,i],digits=d),2d+2))
      println(io,part)
   end; end
   io ≠ stdout ? println(io,"\n\n") : nothing
   io ≠ stdout ? close(io) : nothing
end

function triangleprint(mat;io=stdout,d=4,col=5)
   nams = string.(collect(1:size(mat,1)))
   triangleprint(mat,nams,io=io,d=d,col=col)
end


