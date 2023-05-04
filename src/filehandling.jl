"""
This contains all the filehandling for the westerfit package. It's not particularly well
   written at this point. Very much open for other ideas.
"""
#####INPUTS
function ctrlinit()
   ctrl = Dict("NFOLD" => 0, "S" => 0., "TK" => 8.0, "mcalc" => 8, "vtmax" => 0,
      "Jmax" => 0, "apology" => true, "νmin"=>0.0, "νmax"=>40., "INTthres"=>0.0000, 
      "λlm0"=>0.0001, "RUNmode"=>"ESF", "turducken"=>1, "maxiter"=>60)
   return ctrl
end
function ctrlinp(molnam::String)
   findstr = `grep -n CNTRLS $molnam.inp`
   strln = parse(Int,readchomp(pipeline(findstr,`cut -d : -f1`)))
   findstr = `grep -n '^$' $molnam.inp`
   enln = split(readchomp(pipeline(findstr,`cut -d : -f1`)), "\n")[1]
   len = parse(Int, enln) - strln - 1
   ctrl = ctrlinit()
   file = readdlm(pwd()*"/"*molnam*".inp",'=', skipstart=strln,comments=true,comment_char='#')
   for i in 1:len
      nam = strip(file[i,1])
      val = file[i,2]
      ctrl[nam] = val
   end
   #trl["S"] = Float64(ctrl["S"])
   if ctrl["apology"] == true
      println("Sorry about the name...")
   end
   if iseven(2*ctrl["Jmax"])&&isodd(2*ctrl["S"])
      ctrl["Jmax"] += 0.5
      println("Jmax must be half integer for half interger S. Adding 1/2")
   elseif isodd(2*ctrl["Jmax"])&&iseven(2*ctrl["S"])
      ctrl["Jmax"] += 0.5
      println("Jmax must be integer for interger S. Adding 1/2")
   end
   #println(ctrl)
   return ctrl
end

function secordinit()
   prd = Dict("A" => 1, "B" => 2, "C" => 3, "Dab" => 4, "F" => 12, "ρ" => 13,
      "Vn" => 14, "ϵzz" => 5, "ϵxx" => 6, "ϵyy" => 7, "ϵxz" => 8, "η" => 15,
      "χzz"=> 9, "χxmy"=> 10, "χxz"=> 11)
   return prd
end
function sod2prep(prd::Array{Float64})::Array{Float64}
   out = zeros(15)
   prd[1] += prd[12]*prd[13]^2              # Aeff = A + Fρ²
   out[2] = 0.5*(prd[2] + prd[3])                 #BJ
   out[1] = prd[1] - 0.5*(prd[2] + prd[3])        #BK
   out[3] = 0.25*(prd[2] - prd[3])                #B±
   out[4] = prd[4]                                #Dab
   out[5] = -(prd[5] + prd[6] + prd[7]) / √3.0    #T⁰₀(ϵ)
   out[6] = (2.0*prd[5] - prd[6] - prd[7]) / √6.0 #T²₀(ϵ)
   out[7] = -prd[8]                               #T²₁(ϵ)
   out[8] = (prd[6] - prd[7])*0.5                 #T²₂(ϵ)
   out[9] = prd[9]                                #T²₀(χ)
   out[10] = -√(2.0/3.0)*prd[11]                  #T²₁(χ)
   out[11] = prd[10] / √(6.0)                     #T²₂(χ)
   out[12] = prd[12]*csl                          #F
   out[13] = -2.0*prd[13]*prd[12]*csl             #ρF
   out[14] = prd[14]*0.5*csl                      #V3/2
   out[15] = prd[15]                              #η
   return out
end
function paramrecov(prd::Array{Float64})::Array{Float64}
   out = zeros(15)
   if prd[12] != zero(prd[12])
      out[13] = prd[13]/(-2.0*prd[12])       #ρ
   else
      out[13] = 0.0
   end
   out[1] = prd[1] + prd[2] + out[13]*prd[12]    #A
   out[2] = prd[2] + 2.0*prd[3]                  #B
   out[3] = prd[2] - 2.0*prd[3]                  #C
   out[4] = prd[4]                               #Dab
   out[5] = (prd[5]/√3 - prd[6]/√6)/3.0          #ϵzz
   out[6] = (prd[5]/√3 + 2*prd[6]/√6)/6 + prd[8] #ϵxx
   out[7] = (prd[5]/√3 + 2*prd[6]/√6)/6 - prd[8] #ϵyy
   out[8] = -prd[7]                              #ϵxz
   out[9] = prd[9]                               #χzz
   out[11] = -√(3/2)*prd[10]                     #χxz
   out[10] = √6*prd[11]                          #χxx-χyy
   out[12] = prd[12]/csl                         #F (MHz)
   out[14] = 2.0*prd[14] / csl                   #V3
   out[15] = prd[15]                             #η
   return out
end

function secordinp(molnam::String)
   findstr = `grep -n 2NDORDER $molnam.inp`
   strln = parse(Int,readchomp(pipeline(findstr,`cut -d : -f1`)))
   findstr = `grep -n '^$' $molnam.inp`
   enln = split(readchomp(pipeline(findstr,`cut -d : -f1`)), "\n")[2]
   len = parse(Int, enln) - strln - 1
   #println(len)
   secns = secordinit()
   file = readdlm(pwd()*"/"*molnam*".inp",';', skipstart=strln,comments=true,comment_char='#')
   val = zeros(Float64,15)
   err = zeros(Float64,15)
   for i in 1:len
      nam = strip(file[i,1])
      ind = secns[nam]
      val[ind] = file[i,2]
      err[ind] = file[i,3]
   end
   val = sod2prep(val)
   return val, err
end

parsetuple(T::Type,s::AbstractString) = Tuple(parse.(T, split(s, ',')))
parsetuple(T::Type,s::Int) = s
parsetuple(s::Int) = s
parsetuple(T::Type,s::Float64) = s
function opinp(molnam::String)
   findstr = `grep -n PARAMS $molnam.inp`
   strln = parse(Int,readchomp(pipeline(findstr,`cut -d : -f1`))) + 1
   file = try
      readdlm(pwd()*"/"*molnam*".inp",';', skipstart=strln,comments=true,comment_char='#')
   catch
      zeros(0,0)
   end
   len = size(file,1) #might be wrong index
   if len !=zero(len) #normal behavior for added parameters
   nams = fill("nam",len)
   vals = zeros(Float64,len)
   #vals = Array{Any}(nothing,len)
   errs = zeros(Float64,len)
   #oprs = Array{Any}(nothing,8,len)
   oprs = zeros(Int,8,len)
   stgs = zeros(Int,len)
   col = collect(1:len)
   for i in col
      nams[i] = file[i,1]
      vals[i] = parsetuple(Float64,file[i,2])
      oprs[:,i] = parsetuple.(Int,file[i,4:end-1])
      errs[i] = file[i,3]
      stgs[i] = file[i,end]
   end
   prm0 = col[isequal.(stgs,0)]
   prm1 = col[stgs .!= 0]
   #moved the filtering to the 
   #prm2 = collect(1:len)[isequal.(stgs,2)]
   μset = μproc(vals[prm0], oprs[:,prm0])
   opvls = vals[prm1]#, nams[prm2])
   onams = nams[prm1]#, vals[prm2])
   oerr = errs[prm1]#, errs[prm2])
   oopr = oprs[:,prm1]#, oprs[:,prm2])
   ostg = stgs[prm1]
   return μset, opvls, onams, oerr, oopr, ostg
   else #exception for if none
      println("No higher order operators found. Continuing w/ only H^(2)")
   return nothing, nothing, nothing, nothing, nothing
   end
end

function prmsetter(prm,stg::Array{Int})
   out = copy(prm)
   inds = collect(1:length(prm))[isless.(stg,0)]
   for i in inds
      s = stg[i]
      out[i] *= prm[i+s]
   end
   return out
end


function lineprep(lns,nf,s,vtm)#THIS NEEDS TO BE REWORKED FOR VTM behavior
   #converts the input file into a more code friendly format
   #           1  2  3   4   5  6  7  8   9  10  11   12
   #input  = [ju nu kau kcu mu jl nl kal kcl ml freq unc]
   #           1  2   3   4  5   6
   #output = [ju σu indu jl σl indl]
if nf != zero(nf)
   qunus = lns[:,1:10]
   freqs = lns[:,11]
   uncs = lns[:,12]
   inds = zeros(Int,size(lns,1),6)
   inds[:,1] = Int.(2 .* qunus[:,1])
   inds[:,2] = Int.(mod.(qunus[:,5],nf))
   inds[:,3] = qn2ind.(nf,vtm,qunus[:,5],qunus[:,1],s,qunus[:,2],qunus[:,3],qunus[:,4])
   inds[:,4] = Int.(2 .* qunus[:,6])
   inds[:,5] = Int.(mod.(qunus[:,10],nf))
   inds[:,6] = qn2ind.(nf,vtm,qunus[:,10],qunus[:,6],s,qunus[:,7],qunus[:,8],qunus[:,9])
else
   qunus = lns[:,1:10]
   freqs = lns[:,11]
   uncs = lns[:,12]
   inds = zeros(Int,size(lns,1),6)
   inds[:,1] = Int.(2 .* qunus[:,1])
   inds[:,2] .= 0
   inds[:,3] = qn2ind.(nf,vtm,qunus[:,5],qunus[:,1],s,qunus[:,2],qunus[:,3],qunus[:,4])
   inds[:,4] = Int.(2 .* qunus[:,6])
   inds[:,5] .= 0
   inds[:,6] = qn2ind.(nf,vtm,qunus[:,10],qunus[:,6],s,qunus[:,7],qunus[:,8],qunus[:,9])
end
   return inds, freqs, uncs
end


function pred2lne(sim,s)
"""
Converts the transitions output from westersim into the line format for
   westerfit. Allows for quick test scripting
"""
   if s != zero(s)
      out = zeros(Float64,size(sim)[1],12)
      out[:,1:10] = sim[:,1:10]
      out[:,11] = sim[:,11]
      out[:,12] = fill(0.08,size(sim)[1])
   else#this should be though
      out = zeros(Float64,size(sim)[1],12)
      out[:,1] = sim[:,1]
      out[:,2:5] = sim[:,1:4]
      out[:,6] = sim[:,5]
      out[:,7:10] = sim[:,5:8]
      out[:,11] = sim[:,9]
      out[:,12] = fill(0.08,size(sim)[1])
   end
   return out
end


#####OUTPUTS
function reslin(line,omc,cfrq)
   part  = lpad(line[ 1],6)*";"
   part *= lpad(Int(line[ 2]),4)*";"
   part *= lpad(Int(line[ 3]),4)*";"
   part *= lpad(Int(line[ 4]),4)*";"
   part *= lpad(Int(line[ 5]),4)*";"
   part *= lpad(line[ 6],6)*";"
   part *= lpad(Int(line[ 7]),4)*";"
   part *= lpad(Int(line[ 8]),4)*";"
   part *= lpad(Int(line[ 9]),4)*";"
   part *= lpad(Int(line[10]),4)*";"
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
      part *= @sprintf("%8.4f", frql[2])*","
      part *= @sprintf("%10.4f", frql[3])
   else
      part  = lpad(qunl[1],3)*"/2,"
      part *= lpad(qunl[2],3)*","
      part *= lpad(qunl[3],3)*","
      part *= lpad(qunl[4],3)*","
      part *= lpad(qunl[5],3)*","
      part *= lpad(qunl[7],3)*"/2,"
      part *= lpad(qunl[8],3)*","
      part *= lpad(qunl[9],3)*","
      part *= lpad(qunl[10],3)*","
      part *= lpad(qunl[11],3)*","
      part *= " "*@sprintf("%13.4f", frql[1])*","
      part *= @sprintf("%8.4f", frql[2])*","
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
