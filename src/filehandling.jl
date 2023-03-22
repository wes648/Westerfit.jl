"""
This contains all the filehandling for the westerfit package. It's not particularly well
   written at this point. Very much open for other ideas.
"""
#####INPUTS
function ctrlinit()
   ctrl = Dict("NFOLD" => 0, "S" => 0., "TK" => 8.0, "mcalc" => 8, "vtmax" => 0,
      "Jmax" => 0, "apology" => true, "νmin"=>0.0, "νmax"=>40., "INTthres"=>0.0000, 
      "λ"=>0.0, "RUNmode"=>"ESF")
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
   if ctrl["apology"] == true
      println("Sorry about the name...")
   end
   if iseven(2*ctrl["Jmax"])&&isodd(2*ctrl["S"])
      ctrl["Jmax"] -= 0.5
      println("Jmax must be half integer for half interger S. Subtracting 1/2")
   elseif isodd(2*ctrl["Jmax"])&&iseven(2*ctrl["S"])
      ctrl["Jmax"] -= 0.5
      println("Jmax must be integer for interger S. Subtracting 1/2")
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
   out[1] = prd[1] - out[2]                       #BK
   out[3] = 0.25*(prd[2] - prd[3])                #B±
   out[4] = prd[4]                                #Dab
   out[5] = -(prd[5] + prd[6] + prd[7]) / √3.0    #T⁰₀(ϵ)
   out[6] = (2.0*prd[5] - prd[6] - prd[7]) / √6.0 #T²₀(ϵ)
   out[7] = -prd[8]                               #T²₁(ϵ)
   out[8] = (prd[6] - prd[7])*0.5                 #T²₂(ϵ)
   out[9] = prd[9]                                #T²₀(χ)
   out[10] = -√(2.0/3.0)*prd[11]                  #T²₁(χ)
   out[11] = prd[10] / √6.0                       #T²₂(χ)
   out[12] = prd[12]*csl                          #F
   out[13] = 2.0*prd[13]*prd[12]                  #ρF
   out[14] = prd[14]*0.5*csl                      #V3/2
   out[15] = prd[15]                              #η
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
   errs = zeros(Float64,len)
   oprs = zeros(Int,8,len)
   stgs = zeros(Int,len)
   col = collect(1:len)
   for i in col
      nams[i] = file[i,1]
      vals[i] = file[i,2]
      oprs[:,i] = file[i,4:end-1]
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

function prmsetter(prm::Array{Float64},stg::Array{Int})::Array{Float64}
   out = copy(prm)
   inds = collect(1:length(prm))[isless.(stg,0)]
   for i in inds
      s = stg[i]
      out[i] *= prm[i+s]
   end
   return out
end
#=
function hstager(nams::Array,vals::Array,errs::Array,oprs::Array)
   len = length(val)
   stg = opr[:,end]
   prm1 = collect(1:len)[isequal.(stgs,1)]
   prm2 = collect(1:len)[isequal.(stgs,2)]
   prm3 = collect(1:len)[isequal.(stgs,3)]
   onam = (nothing, nothing, nothing)
   oval = (nothing, nothing, nothing)
   oerr = (nothing, nothing, nothing)
   oopr = (nothing, nothing, nothing)
   if length(prm1) != 0
      onam[1] = nams[prm1]
      oval[1] = vals[prm1]
      oerr[1] = errs[prm1]
      oopr[1] = oprs[prm1]
   end
   if length(prm2) != 0
      onam[2] = nams[prm2]
      oval[2] = vals[prm2]
      oerr[2] = errs[prm2]
      oopr[2] = oprs[prm2]
   end
   if length(prm3) != 0
      onam[3] = nams[prm3]
      oval[3] = vals[prm3]
      oerr[3] = errs[prm3]
      oopr[3] = oprs[prm3]
   end
   return onam, oval, oerr, oopr
end
function hstager(nams::Nothing,vals::Nothing,errs::Nothing,oprs::Nothing)
   out = (nothing, nothing, nothing)
   return out, out, out, out
end

function hstager_sim(vals::Array,oprs::Array)
   len = length(vals)
   stg = oprs[end,:]
   prm1 = collect(1:len)[isequal.(stg,1)]
   prm2 = collect(1:len)[isequal.(stg,2)]
   prm3 = collect(1:len)[isequal.(stg,3)]
   v1 = nothing
   p1 = nothing
   v2 = nothing
   p2 = nothing
   v3 = nothing
   p3 = nothing
   if length(prm1) != 0
      v1 = vals[prm1]
      p1 = oprs[:,prm1]
   end
   if length(prm2) != 0
      v2 = vals[prm2]
      p2 = oprs[:,prm2]
   end
   if length(prm3) != 0
      v3 = vals[prm3]
      p3 = oprs[:,prm3]
   end
   oval = (v1, v2, v3)
   oopr = (p1, p2, p3)
   return oval, oopr
end
function hstager_sim(vals::Nothing,oprs::Nothing)
   out = (nothing, nothing, nothing)
   return out, out, out, out
end
=#
function paraminit()
   prd = Dict("A" => 1, "B" => 2, "C" => 3, "Dab" => 4, "F" => 5, "ρ" => 6,
      "V3" => 7, "ϵzz" => 8, "ϵxx" => 9, "ϵyy" => 10, "ϵxz" => 11, "η" => 12,
      "χzz" => 13, "χxmy" => 14, "χxz" => 15, "ΔN" => 16, "ΔNK" => 17, "ΔK" => 18,
      "δN" => 19, "δK" => 20, "Fm" => 21, "V6" => 22, "V3m" => 23, "ρm" => 24,
      "ρ3" => 25, "FN" => 26, "FK" => 27, "Fbc" => 28, "Fab" => 29, "V3N" => 30,
      "V3K" => 31, "V3ab" => 32, "V3bc" => 33, "ρN" => 34, "ρK" => 35, "ρab" => 36,
      "ρbN" => 37, "ΔsN" => 38, "ΔsNK" => 39, "ΔsKN" => 40, "ΔsK" => 41, "δsN" => 42,
      "δsK" => 43, "ΦJ" => 44, "ΦJK" => 45, "ΦKJ" => 46, "ΦK" => 47, "ϕJ" => 48,
      "ϕJK" => 49, "ϕK" => 50, "μa" => 51, "μb" => 52, "μc" => 53)
   return prd
end
function rotfix!(prm)
   BJ = 0.5*(prm[2]+prm[3])
   BK = prm[1] - BJ
   Bp = 0.5*(prm[2]-prm[3])
   prm[1] = BK
   prm[2] = BJ
   prm[3] = Bp
   return prm
end
function spifix!(prm)
   ao = -(prm[8] + prm[9] + prm[10])/√3.0
   a = -(2.0*prm[8] - prm[9] - prm[10])/√6.0
   d = -prm[11]*0.5
   b = (prm[9] - prm[11])*0.5
   χ2 = √(1.0/6.0)*prm[14]
   χ1 = -√(2.0/3.0)*prm[15]
   prm[8] = ao
   prm[9] = a
   prm[10] = b
   prm[11] = d
   prm[14] = χ2
   prm[15] = χ1
   return prm
end
function paraminp(molnam::String)
   findstr = `grep -n PARAMS $molnam.inp`
   strln = parse(Int,readchomp(pipeline(findstr,`cut -d : -f1`)))
   inds = paraminit()
   prms = zeros(53)
   errs = zeros(53)
   file = readdlm(pwd()*"/"*molnam*".inp",',', skipstart=strln)
   for i in 1:size(file,1)
      ind = inds[file[i,1]]
      val = file[i,2]
      err = file[i,3]
      prms[ind] = val
      errs[ind] = err
   end
   prms = rotfix!(prms)
   prms[6] *= prms[5]
   prms = spifix!(prms)
   μs = [prms[end-2:end] errs[end-2:end]]
   prms = prms[1:end-3]
   errs = errs[1:end-3]
   return prms, errs, μs
end


function lineprep(lns,nf,s,mcalc)#THIS NEEDS TO BE REWORKED
   #           1  2  3   4   5  6  7  8  9   10 11 12  13   14
   #input  = [ju nu kau kcu mu σu jl nl kal kcl ml σl freq unc]
   #converts the input file into a more code friendly format
   #           1  2  3   4   5  6  7  8   9  10  11   12
   #input  = [ju nu kau kcu mu jl nl kal kcl ml freq unc]
   #           1  2   3   4  5   6
   #output = [ju σu indu jl σl indl]
   qunus = lns[:,1:10]
   freqs = lns[:,11]
   uncs = lns[:,12]
   inds = zeros(Int,size(lns,1),6)
   inds[:,1] = Int.(2 .* qunus[:,1])
   inds[:,2] = Int.(mod.(qunus[:,5],nf))
   inds[:,3] = qn2ind.(nf,mcalc,qunus[:,5],qunus[:,1],s,qunus[:,2],qunus[:,3],qunus[:,4])
   inds[:,4] = Int.(2 .* qunus[:,6])
   inds[:,5] = Int.(mod.(qunus[:,10],nf))
   inds[:,6] = qn2ind.(nf,mcalc,qunus[:,10],qunus[:,6],s,qunus[:,7],qunus[:,8],qunus[:,9])
   #inds = vcat(inds[:,1:2], inds[:,3:4])
   return inds, freqs, uncs
end


#####OUTPUTS
function englin(s,eng,qunl)
   if s==zero(s)
      part = lpad(qunl[2],4)*","
   else
      part = lpad(qunl[1],4)*"/2,"
      part *= lpad(qunl[2],4)*","
   end
   part *= lpad(qunl[3],4)*","
   part *= lpad(qunl[4],4)*","
   part *= lpad(qunl[5],4)*","
   part *= lpad(qunl[6],4)*","
   part *= " "*lpad(@sprintf("%0.10f", eng), 16)
   return part
end

function EngWriter(molnam,ctrl,energies,qunus)
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
