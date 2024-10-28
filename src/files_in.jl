"""
This contains all the input file handling for the westerfit package. It's not 
   particularly well written at this point. Very much open for other ideas.
"""
#####INPUTS
function ctrlinit()
   ctrl = Dict("NFOLD" => 0, "S" => 0., "TK" => 8.0, "mcalc" => 8, "vtmax" => 0,
      "Jmax" => 0, "apology" => true, "νmin"=>0.0, "νmax"=>40., "INTthres"=>0.00001, 
      "λlm0"=>0.0001, "RUNmode"=>"ESF", "turducken"=>1, "maxiter"=>60, "overwrite"=>true,
      "assign"=>"ram36", "REJECT"=>1.0e+1, "Irrep"=>"Ir", "goal"=>1.0, "mmax"=>6, "stages"=>1)
   return ctrl
end
function blockfind(molnam::String,blknam::String)
   out = zeros(Int,2)
   count = 0
   inblock = false
   io = open("$molnam.inp")
   for i in eachline(io)
      count += 1
      if contains(i,blknam)
         out[1] = count
         inblock = true
      elseif isempty( filter(x -> !isspace(x), i) )&& inblock==true
         out[2] = count
         break
      elseif eof(io)
         out[2] = count
         break
      end
   end
   close(io)
   #@show blknam
   #@show out
   return out
end

function ctrlinp(molnam::String)
   blk = blockfind(molnam, "%CNTRLS")
   len = blk[2] - blk[1] - 1
   ctrl = ctrlinit()
   file = readdlm(pwd()*"/"*molnam*".inp",'=', skipstart=blk[1],comments=true,comment_char='#')
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
   if ctrl["NFOLD"]==0
      ctrl["mcalc"] = 0
      ctrl["vtmax"] = 0
   end
   ctrl["assign"] = strip(ctrl["assign"])
   ctrl["Irrep"] = String(strip(ctrl["Irrep"]))
   #println(ctrl)
   return ctrl
end

function secordinit_old()
   prd = Dict("A" => 1, "B" => 2, "C" => 3, "Dab" => 4, "Dxz" => 4, "F" => 12, "ρ" => 13,
      "Vn" => 14, "ϵzz" => 5, "ϵxx" => 6, "ϵyy" => 8, "ϵxz" => 7, "η" => 15,
      "χzz"=> 9, "χxmy"=> 11, "χxz"=> 10)
   return prd
end
function secordinit_full()
   prd = Dict("A" => 1, "B" => 2, "C" => 3, "Dab" => 4, "Dxz" => 4, 
      "ϵzz" => 5, "ϵxx" => 6, "ϵyy" => 7, "ϵzx" => 8, "ϵxz" => 9,
      "Czz" => 5, "Cxx" => 6, "Cyy" => 7, "Czx" => 8, "Cxz" => 9,
      "χzz"=> 10, "χxz"=> 11, "χxmy"=> 12, "χxx-χyy"=>12,
      "α"=>10, "δ"=>11, "β"=>12,
      "F" => 13, "ρz" => 14, "ρ" => 14, "ρx" =>15, "Vn" => 16, "V3" =>16,
      "ηz" => 17, "η" => 17, "ηx" => 18)
   return prd
end

function irrepswap(irrep::String)::Array{Int}
   if irrep=="Il"
      perm = [1;3;2]
   elseif irrep=="IIr"
      perm = [2;3;1]
   elseif irrep=="IIl"
      perm = [2;1;3]
   elseif irrep=="IIIr"
      perm = [3;1;2]
   elseif irrep=="IIIl"
      perm = [3;2;1]
   else #Ir is default
      perm = [1;2;3]
   end
   return perm
end

function sod2prep_old(prd::Array{Float64})::Array{Float64}
   out = zeros(15)
   tempa = prd[1] + csl*prd[12]*prd[13]^2         #Aeff = A + Fρ²
   out[2] = 0.5*(prd[2] + prd[3])                 #BN
   out[1] = tempa - 0.5*(prd[2] + prd[3])         #BK
   out[3] = 0.25*(prd[2] - prd[3])                #B±
   out[4] = prd[4]                                #Dab
   out[5] = -(prd[5] + prd[6] + prd[8]) / √3      #T⁰₀(ϵ)
   out[6] = (2*prd[5] - prd[6] - prd[8]) / √6     #T²₀(ϵ) 
   out[7] = -prd[7]                               #T²₁(ϵ)
   out[8] = (prd[6] - prd[8])*0.5                 #T²₂(ϵ)
   out[9] = prd[9]                                #T²₀(χ)
   out[10] = -√(2.0/3.0)*prd[10]                  #T²₁(χ)
   out[11] = prd[11] / √(6.0)                     #T²₂(χ)
   out[12] = prd[12]*csl                          #F
   out[13] = -2.0*prd[13]*prd[12]*csl             #ρF
   out[14] = prd[14]*0.5*csl                      #V3/2
   out[15] = prd[15]                              #η
   return out
end
function sod2prep_full(prd::Array{Float64})::Array{Float64}
   out = zeros(18)
   tempa = prd[1] + csl*prd[13]*prd[14]^2         #Aeff = A + Fρx²
   tempb = prd[2] + csl*prd[13]*prd[15]^2         #Beff = B + Fρz²
   out[1] = tempa - 0.5*(tempb + prd[3])          #BK
   out[2] = 0.5*(tempb + prd[3])                  #BN
   out[3] = 0.25*(tempb - prd[3])                 #B±
   out[4] = prd[4]                                #Dab

   out[5] = -(prd[5] + prd[6] + prd[7]) / √3      #T⁰₀(ϵ)
   out[6] = 0.5*(prd[8] - prd[9])                 #T¹₁(ϵ)
   out[7] = (2*prd[5] - prd[6] - prd[7]) / √6     #T²₀(ϵ)
   out[8] = -0.5*(prd[8] + prd[9])                #T²₁(ϵ)
   out[9] = (prd[6] - prd[7])*0.5                 #T²₂(ϵ)
   
   out[10] = prd[10]                              #T²₀(χ)
   out[11] = -√(2.0/3.0)*prd[11]                  #T²₁(χ)
   out[12] = prd[12] / √(6.0)                     #T²₂(χ)
   #the factor of 0.5 on the x terms is from 2lₓ = l_+ + l_-
   out[13] = prd[13]*csl                          #F
   out[14] = -2.0*prd[13]*prd[14]*csl             #ρzF
   out[15] = -prd[13]*prd[15]*csl                 #ρxF
   out[16] = prd[16]*0.5*csl                      #Vn/2
   out[17] = prd[17]                              #ηz
   out[18] = 0.5*prd[18]                          #ηx
   return out
end

function epszxcheck!(prd::Array{Float64},scl::Array{Float64})
   if prd[6] ≈ prd[8]
      prd[6] = 0.0
      scl[6] = 0.0
      prd[8] *= 2.0
   end
   return prd, scl
end

function secordinp(molnam::String,irp::String)
   blk = blockfind(molnam, "%2NDORDER")
   len = blk[2] - blk[1] - 1
   #println(len)
   secns = secordinit_full()
   file = readdlm(pwd()*"/"*molnam*".inp",';', skipstart=blk[1])#,comments=true,comment_char='#')
   #println(file)
   val = zeros(Float64,19)
   err = zeros(Float64,19)
   for i in 1:len
      nam = strip(file[i,1])
      #ind = secns[nam]
      ind = get(secns, nam, 19)
      val[ind] = file[i,2]
      err[ind] = file[i,3]
   end
   val[1:3] = val[irrepswap(irp)]
   pop!(val)
   pop!(err)
   val = sod2prep_full(val)
   val, err = epszxcheck!(val,err)
   #@show val
   return val, err
end

parsetuple(T::Type,s::AbstractString) = Tuple(parse.(T, split(s, ',')))
parsetuple(T::Type,s::Int) = s
parsetuple(s::Int) = s
parsetuple(T::Type,s::Float64) = s
function opreadconv(x)::Int
   if typeof(x) == SubString{String}
      return 0
   else
      return x
   end
end

function opinp(molnam::String,nfold)
   #findstr = `grep -n PARAMS $molnam.inp`
   #strln = parse(Int,readchomp(pipeline(findstr,`cut -d : -f1`))) + 1
   blk = blockfind(molnam, "%PARAMS")
   blk[1] += 1
   len = blk[2] - blk[1] + 1
   file = try
      readdlm(pwd()*"/"*molnam*".inp",';', skipstart=blk[1],comments=true,comment_char='#')
   catch
      zeros(0,0)
   end
   len = size(file,1) #might be wrong index
   if size(file,2)==12
   @warn "You are using a deprecated input format! Check the github for the new structure"
   end
   if len !=zero(len) #normal behavior for added parameters
   nams = fill("nam",len)
   vals = zeros(Float64,len)
   #vals = Array{Any}(nothing,len)
   errs = zeros(Float64,len)
   #oprs = Array{Any}(nothing,8,len)
   oprs = zeros(Int,9,len)
   stgs = zeros(Int,len)
   col = collect(1:len)
   for i in col
      nams[i] = file[i,1]
      vals[i] = file[i,2] #parsetuple(Float64,file[i,2])
      oprs[:,i] = opreadconv.(file[i,4:end-1]) #parsetuple.(Int,file[i,4:end-1])
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
   if (nfold ≠ 0)&&(isodd(nfold))
      oopr[8,:] = floor.(Int,oopr[8,:]./nfold)
      oopr[9,:] = floor.(Int,oopr[9,:]./nfold)
   elseif (nfold ≠ 0)&&(isodd(nfold))
      oopr[8,:] = floor.(Int,2 .*oopr[8,:]./nfold) #dividing by 0.5nfold is the same
      oopr[9,:] = floor.(Int,2 .*oopr[9,:]./nfold) # as multipling by 2/nfold
   end
   ostg = stgs[prm1]
   return μset, opvls, onams, oerr, oopr, ostg
   else #exception for if none
      println("No higher order operators found. Continuing w/ only H^(2)")
   return nothing, nothing, nothing, nothing, nothing
   end
end

function prmsetter(prm,stg::Array{Int})
   #This function runs inside the TSRCALC functions
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
if length(lns[1,:]) != 12
    println("Wrong number of columns in .lne file")
end
if nf≠zero(nf) && isodd(nf)
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
elseif nf≠zero(nf) && iseven(nf)
   qunus = lns[:,1:10]
   freqs = lns[:,11]
   uncs = lns[:,12]
   inds = zeros(Int,size(lns,1),6)
   inds[:,1] = Int.(2 .* qunus[:,1])
   inds[:,2] = Int.(mod.(2 .*qunus[:,5],nf))
   inds[:,3] = qn2ind.(nf,vtm,qunus[:,5],qunus[:,1],s,qunus[:,2],qunus[:,3],qunus[:,4])
   inds[:,4] = Int.(2 .* qunus[:,6])
   inds[:,5] = Int.(mod.(2 .*qunus[:,10],nf))
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

function qnconv(lns,nf,s,vtm)#THIS NEEDS TO BE REWORKED FOR VTM behavior
   #converts the simulator quantum numbers into the derivative format
   #           1  2  3   4   5  6  7  8   9  10 
   #input  = [ju nu kau kcu mu jl nl kal kcl ml]
   #           1  2   3   4  5   6
   #output = [ju σu indu jl σl indl]
if length(lns[1,:]) != 10
    println("Wrong number of columns in qns")
end
if nf ≠ zero(nf)
   qunus = lns[:,1:10]
   inds = zeros(Int,size(lns,1),6)
   inds[:,1] = Int.(2 .* qunus[:,1])
   inds[:,2] = Int.(mod.(qunus[:,5],nf))
   inds[:,3] = qn2ind.(nf,vtm,qunus[:,5],qunus[:,1],s,qunus[:,2],qunus[:,3],qunus[:,4])
   inds[:,4] = Int.(2 .* qunus[:,6])
   inds[:,5] = Int.(mod.(qunus[:,10],nf))
   inds[:,6] = qn2ind.(nf,vtm,qunus[:,10],qunus[:,6],s,qunus[:,7],qunus[:,8],qunus[:,9])
else
   qunus = lns[:,1:10]
   inds = zeros(Int,size(lns,1),6)
   inds[:,1] = Int.(2 .* qunus[:,1])
   inds[:,2] .= 0
   inds[:,3] = qn2ind.(nf,vtm,qunus[:,5],qunus[:,1],s,qunus[:,2],qunus[:,3],qunus[:,4])
   inds[:,4] = Int.(2 .* qunus[:,6])
   inds[:,5] .= 0
   inds[:,6] = qn2ind.(nf,vtm,qunus[:,10],qunus[:,6],s,qunus[:,7],qunus[:,8],qunus[:,9])
end
   return inds
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


