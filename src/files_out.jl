"""
This handles the output files for the westerfit package. Sophie made most of this!
"""

function paramrecov(prd::Array{Float64})::Array{Float64}
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
   oprd[1:3] = oprd[irrepswap(irrep)]
   ounc[1:3] = ounc[irrepswap(irrep)]
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
   eng = energies[:,1]
   qns = qunus[:,:,1]
   for sc in 2:σcount(ctrl["NFOLD"])
      eng = vcat(eng,energies[:,sc])
      qns = vcat(qns,qunus[:,:,sc])
   end
   len = size(eng,1)
   out = fill("0",len)
   for i in 1:len
      energy = eng[i]/c
      #0.10f is what BELGI uses, 0.6f is for spcat
      out[i] = englin(ctrl["S"],energy,qns[i,:])
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
   for sc in 2:σcount(ctrl["NFOLD"])
      eng = vcat(eng,energies[:,sc])
      qns = vcat(qns,qunus[:,:,sc])
      psz = vcat(psz,pasz[:,sc])
   end
   len = size(eng,1)
   out = fill("0",len)
   for i in 1:len
      energy = eng[i]/c
      #0.10f is what BELGI uses, 0.6f is for spcat
      out[i] = englin(ctrl["S"],energy,qns[i,:],pasz[i])
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
      part  = lpad(qunl[1]*0.5,4)*","
      part *= lpad(qunl[2],3)*","
      part *= lpad(qunl[3],3)*","
      part *= lpad(qunl[4],3)*","
      part *= lpad(qunl[5],3)*","
      part *= lpad(qunl[7]*0.5,4)*","
      part *= lpad(qunl[8],3)*","
      part *= lpad(qunl[9],3)*","
      part *= lpad(qunl[10],3)*","
      part *= lpad(qunl[11],3)*","
      part *= " "*@sprintf("%13.4f", frql[1])*","
      part *= @sprintf("%10.4f", frql[4])*","
      part *= @sprintf("%12.6f", frql[2])*","
      part *= @sprintf("%10.4f", frql[3])
   end
   return part
end
function TraWriterSPCAT(molnam,freqs, qunus) #emulates the cat file structure of SPCAT
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
function TraWriter(molnam,s, freqs, qunus)
   p = sortperm(freqs[:,1])
   freqs = freqs[p,:]
   qunus = qunus[p,:]
   out = fill("0",size(freqs,1))
   for i in 1:size(freqs,1)
      out[i] = linestrng(s, freqs[i,:], qunus[i,:])
   end
   io = open("$molnam.cat", "w") do io
      for i in out
         println(io, i)
      end
   end
   println("Transitions written to $molnam.cat!")
end


function TraWriterold(molnam,freqs, qunus)
   c = 29979.2458
   p = sortperm(freqs[:,1])
   freqs = freqs[p,:]
   qunus = qunus[p,:]
   out = fill("0",size(freqs,1))
   counter = 0
   for i in 1:size(freqs,1)
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


function inpwriter(molnam::String, values, scales)

   iterativecopier(molnam)

   strlnctrl,strln2nd,strlnhigh,file = findstrinput_new(molnam)

   #strlnctrl = first(findall(isequal("%CNTRLS"),file[:,1]))
   #strln2nd = first(findall(isequal("%2NDORDER"),file[:,1]))
   #strlnhigh = first(findall(isequal("%PARAMS NᵃSᵇNzᶜSzᵈ(N₊ᵉS₊ᶠ + S₋ᶠN₋ᵉ)Pₐᵍcos(hα)sin(jα)Ny^(1-δ(0,j))"),file[:,1]))

  # numvals = strlnhigh-strln2nd-1
   secvalues = values[1:18] #new values of second order ops
   #println(secvalues)

   highstg= file[strlnhigh+2:end,end] #stages of higher order ops
   ohighval = file[strlnhigh+2:end,2] #old vals of higher order ops
   uval = ohighval[highstg .== 0.0] #values only of dipole moments
   newhighval = values[19:end] #new values of higher order ops

   highervalues = vcat(uval,newhighval) #all higher order ops, incl dipoles


   controls = file[strlnctrl:strln2nd-1, 1]
#   secnam = file[strln2nd+1:strlnhigh-1,1]
#   secscale = file[strln2nd+1:strlnhigh-1,3]

   secnam = ["A","B","C","Dab","ϵzz","ϵxx","ϵyy","ϵzx","ϵxz","χzz","χxz","χxmy","F","ρz","ρx","Vn","ηz","ηx"]

   highnam = file[strlnhigh+2:end,1]
   highscale = file[strlnhigh+2:end,3]
   higha= file[strlnhigh+2:end,4]
   highb= file[strlnhigh+2:end,5]
   highc= file[strlnhigh+2:end,6]
   highd= file[strlnhigh+2:end,7]
   highe= file[strlnhigh+2:end,8]
   highf= file[strlnhigh+2:end,9]
   highg= file[strlnhigh+2:end,10]
   highh= file[strlnhigh+2:end,11]
   highj= file[strlnhigh+2:end,12]
   
   secondord = fill("0",18)
   higherord = fill("0",length(file[strlnhigh+2:end,1]))


   for i in 1:length(secvalues)
      ln = 30 - length(secnam[i])
      secondord[i] = string(secnam[i],"; ", lpad(secvalues[i],ln),";",lpad(scales[i],6))
   end

   for i in 1:length(higherord)
      ln = 30 - length(highnam[i])
      higherord[i] = string(highnam[i],"; ",lpad(highervalues[i],ln),";",
         lpad(highscale[i],6),";",lpad(higha[i],4),";",lpad(highb[i],4),";",
         lpad(highc[i],4),";",lpad(highd[i],4),";",lpad(highe[i],4),";",
         lpad(highf[i],4),";",lpad(highg[i],4),";",lpad(highh[i],4),";",lpad(highj[i],4),";",lpad(highstg[i],4))
   end

   time = now()
   io = open("$molnam.inp", "w") do io
      println(io,molnam, "   @   ", time)
      for i in 1:length(controls)
         println(io, controls[i])
      end
      println(io, "")
      println(io,"%2NDORDER")
      for i in 1:length(secondord)
         println(io, secondord[i])
      end
      println(io,"")
      println(io,"%PARAMS NᵃSᵇNzᶜSzᵈ(N₊ᵉS₊ᶠ + S₋ᶠN₋ᵉ)Pₐᵍcos(hα)sin(jα)Ny^(1-δ(0,j))")
      println(io,"%Op;                         Val;   Scl;   a;   b;   c;   d;   e;   f;   g;   h;   j;  stg")
      for i in 1:length(higherord)
         println(io, higherord[i])
      end
   end
end

"""END CODE FOR INPUT FILE WRITER"""

function findstrinput(molnam)
   #findctrl = `grep -n CNTRLS $molnam.inp`
   #strlnctrl = parse(Int,readchomp(pipeline(findctrl,`cut -d : -f1`)))
   #find2nd = `grep -n 2NDORDER $molnam.inp`
   #strln2nd = parse(Int,readchomp(pipeline(find2nd,`cut -d : -f1`)))
   #findhigh = `grep -n PARAMS $molnam.inp`
   #strlnhigh = parse(Int,readchomp(pipeline(findhigh,`cut -d : -f1`)))
   strlnctrl = blockfind(molnam,"CNTRLS")[1]
   strln2nd = blockfind(molnam,"2NDORDER")[1]
   strlnhigh = blockfind(molnam,"PARAMS")[1]
   file = readdlm(pwd()*"/"*molnam*".inp",';',comments=true,comment_char='#')
   return strlnctrl,strln2nd,strlnhigh,file
end

function findstrinput_new(molnam)
   file = readdlm(pwd()*"/"*molnam*".inp",';',comments=true,comment_char='#')
   strlnctrl = findfirst(isequal("%CNTRLS"),file)[1]
   strln2nd = findfirst(isequal("%2NDORDER"),file)[1]
   strlnhigh = findfirst(isequal("%PARAMS NᵃSᵇNzᶜSzᵈ(N₊ᵉS₊ᶠ + S₋ᶠN₋ᵉ)Pₐᵍcos(hα)sin(jα)Ny^(1-δ(0,j))"),file)[1]
   return strlnctrl,strln2nd,strlnhigh,file
end


function outputinit(molnam,params,scls,linelength,ctrl)

  strlnctrl,strln2nd,strlnhigh,file = findstrinput(molnam)

   indexcon = first(findall(isequal("%CNTRLS"),file[:,1]))
   index2nd = first(findall(isequal("%2NDORDER"),file[:,1]))

   controls = file[indexcon+1:index2nd-2, 1]
   
   secvalues = params[1:18]
   highervalues = params[19:end]
   secscale = scls[1:18]
   highscale = scls[19:end]
   highstg= file[strlnhigh:end,end]


   secnam = [" BK", " BN", " B⨦", " Dab", " T⁰₀(ϵ)", " T¹₁(ϵ)"," T²₀(ϵ)",
             " T²₁(ϵ)"," T²₂(ϵ)", " T²₀(χ)"," T²₁(χ)"," T²₂(χ)", 
             " F", " -2ρzF", " -ρxF", " V3/2", " ηz", " ηx"]
   highnamall = file[strlnhigh:end,1]
   highnam = highnamall[highstg .!= 0.0]
   fullnam = vcat(secnam, highnam)
   
   secondord = fill("0",18)
   higherord = fill("0",length(highnam))
   
   for i in 1:length(secondord)
      ln = 30 - length(secnam[i])
      secondord[i] = string(secnam[i],"; ", lpad(secvalues[i],ln),";",lpad(secscale[i],6))
   end
   
   for i in 1:length(higherord)
      ln = 30 - length(highnam[i])
      higherord[i] = string(" ",highnam[i],"; ",lpad(highervalues[i],ln),";",lpad(highscale[i],6))
   end
   
   time = now()
   
   io = open("$molnam.out","w")
      println(io,molnam,"   @   ",time)
      println(io,"")
      println(io,"Control Parameters")
      for (key,value) in ctrl
         println(io,key," = ", value)
      end
      println(io,"")
      println(io,"Initial Parameters")
      for i in 1:length(secvalues)
         println(io,secondord[i])
      end
      println(io,"")
      for i in 1:length(higherord)
         println(io,higherord[i])
      end
      println(io,"")
      println(io,"Number of lines = $linelength")
      println(io,"")
   close(io)
end


function iterationwriter(molnam,paramarray,rms,counter,λlm,βf,perm)

   strlnctrl,strln2nd,strlnhigh,file = findstrinput(molnam)

   #counter = parse(Int,scounter)
   scounter = lpad(counter,3)
   srms = (@sprintf("%0.4f", rms))
   slλ = (@sprintf("%0.4f", log10(λlm)))


   prd = paramarray[:,counter+1]
   prdold = paramarray[:,counter]
   highstg= file[strlnhigh:end,end]

   
   highervalues = prd[19:end]
   
   secnam = [" BK", " BN", " B⨦", " Dab", " T⁰₀(ϵ)", " T¹₁(ϵ)"," T²₀(ϵ)",
             " T²₁(ϵ)"," T²₂(ϵ)", " T²₀(χ)"," T²₁(χ)"," T²₂(χ)", 
             " F", " -2ρzF", " -ρxF", " V3/2", " ηz", " ηx"]
   highnamall = file[strlnhigh:end,1]
   highnam = highnamall[highstg .!= 0.0]
   fullnam = vcat(secnam, highnam)
   
   change = zeros(length(prd))
   change[perm] .+= βf

   schange = fill("Fixed",length(change))
   for i in 1:length(change)
      if change[i] != 0
        schange[i] = @sprintf("%0.4f",change[i])
      else
      end
   end
   
   percent = zeros(length(change))
   for i in 1:length(change)
      if prdold[i] !=0   
         percent[i] = abs((change[i]/prdold[i])*100)
      else
         percent[i] = 0.0
      end
   end
   maxloc = argmax(percent)
   spercent = @sprintf("%0.4f", percent[maxloc])
   maxname = fullnam[maxloc]
   

   prdprint = fill("0",length(fullnam))
   for i in 1:length(prdprint)
      ln = 30 - length(fullnam[i])
      prdprint[i] = string(" ",fullnam[i],"; ", lpad(prd[i],ln), "; ",lpad(schange[i],10))
   end
      
   io = open("$molnam.out", "a")
      println(io, "After $scounter Iterations:")
      println(io, "")
      println(io, "RMS = $srms MHz, log₁₀(λ) = $slλ")
      println(io, "")
      println(io, "                        Parameter     Change")
      for i in 1:length(prdprint)
         println(io, prdprint[i])
      end
      println(io,"")
      println(io,"Max Change; $maxname, $spercent%")
      println(io,"")
      println(io,"-------------------------------------")
      println(io,"")
   close(io)
end

function pdiff(params,puncs)
   pdiff = sfdiv(puncs, params)
return pdiff
end

function outputfinal(molnam,ctrl,frms,counter,slλ,puncs,params,endpoint)

   strlnctrl,strln2nd,strlnhigh,file = findstrinput(molnam)

   srms = (@sprintf("%0.4f", frms))
   scounter = lpad(counter,3)


   secnam = ["A","B","C","Dab","ϵzz","ϵxx","ϵzx","ϵxz","ϵyy","χzz","χxz",
             "χxx-χyy","F","ρz", "ρx", "V3","ηz","ηx"]
   highnamall = file[strlnhigh:end,1]
   highstg= file[strlnhigh:end,end]
   highnam = highnamall[highstg .!= 0.0]
   fullnam = vcat(secnam, highnam)

   io = open("$molnam.out", "a")
   if counter==0
      #This is a super sloppy bug fix for a bug when zero iterations occur
      #println(io,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH")
      #println(io,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH")
      println(io,"hi!! this is to fix a bug in which the first 111 characters appended to a file vanish I don't have a better fix")
   end


      if endpoint == "converge"
         println(io," A miracle has come to pass. The fit has converged.")
      elseif endpoint == "RMS"
         println(io," The RMS has stopped decreasing. Hopefully it is low.")
      elseif endpoint == "step size"
         println(io," It would appear step size has converged.")
      elseif endpoint == "LMthresh"
         println(io," λlm exceeded threshold.")
         println(io," If you were using the turducken, try again without it.")
      elseif endpoint == "iter"
         println(io," Alas, the iteration count has exceeded the limit.")
      else
         println(" The output writer seems to be having issues.")
      end

      println(io,"")

      println(io,"After $scounter iterations:")
      println(io,"")
      println(io,"Final RMS = $srms MHz")
      println(io,"Final log₁₀(λ) = $slλ")
      println(io,"\n")

      unformat = fill("0",length(params))
      println(io,"                         Parameter                   Uncertainty             Percent Difference")
      for i in 1:length(params)
         pdiffi = sfdiv(puncs[i],params[i])*100
         ln = 30 - length(fullnam[i])
         ln2 = 30 - length(params[i])
         ln3 = 30 - length(pdiffi)
         unformat[i] = string(" ",fullnam[i], "; ", lpad(params[i], ln), "; ", lpad(puncs[i], ln2), "; ", lpad(pdiffi, ln3))
      end
      for i in 1:length(unformat)
         println(io,unformat[i])
      end
      println(io,"\n")

      try
         finalunc = uncrformattersci(params,puncs)
         formatted = fill("0",2,length(fullnam))
         for i in 1:length(fullnam)
            ln = 30 - length(fullnam[i])
            formatted[i] = string(fullnam[i], "; ", lpad(finalunc[i],ln))
         end
         println(io,"Journal Formatted Values")
         for i in 1:length(fullnam)
            println(io," ",formatted[i])
         end
      catch
         println(io," Yikes! The uncertainty formatter failed")
         println(io," There are likely some ill-fit parameters")
         println(io," You'll have to manually typeset the uncertainties\n")
      end

   println(io,"\n")
   if ctrl["apology"]==true
      println(io," Again, sorry about the name...")
   end
   println(io,"\n")
   close(io)

   println("Output written to $molnam.out!")
end
