
function ctrlinit()
   ctrl = Dict("NFOLD" => 0, "S" => 0., "TK" => 8.0, "mcalc" => 8, "vtmax" => 0,
      "Jmax" => 0, "apology" => true, "νmin"=>0.0, "νmax"=>40., "INTthres"=>0.00001, 
      "λlm0"=>0.0001, "RUNmode"=>"ESF", "turducken"=>1, "maxiter"=>60, "overwrite"=>true,
      "assign"=>"ram36", "REJECT"=>1.0e+1, "Irrep"=>"Ir", "goal"=>1.0, "mmax"=>6, "stages"=>1,
      "ctbk"=>[0;0],"sobk"=>[0;0],"opbk"=>[0;0])
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
   ctrl = ctrlinit()
   blk = blockfind(molnam, "%CNTRLS")
   len = blk[2] - blk[1] - 1
   ctrl["ctbk"] = blk
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

function secordinit_lim()::Dict{String,Int}
   prd = Dict("A" => 1, "B" => 2, "C" => 3, "Dab" => 4, "Dxz" => 4, 
      "ϵzz" => 5, "ϵxx" => 6, "ϵyy" => 7, "ϵzx" => 8, "ϵxz" => 8,
      "Czz" => 5, "Cxx" => 6, "Cyy" => 7, "Czx" => 8, "Cxz" => 8,
      "χzz"=> 9, "χxz"=> 10, "χxmy"=> 11, "χxx-χyy"=>11,
      "α"=>9, "δ"=>10, "β"=>11)#,
      #"F" => 13, "ρz" => 14, "ρ" => 14, "ρx" =>15, "Vn" => 16, "V3" =>16,
      #"ηz" => 17, "η" => 17, "ηx" => 18)
   return prd
end
function sod2prep_lim(prd::Array{Float64})::Array{Float64}
   out = zeros(12)
   ##tempa = prd[1] + csl*prd[13]*prd[14]^2         #Aeff = A + Fρz²
   ##tempb = prd[2] + csl*prd[13]*prd[15]^2         #Beff = B + Fρx²
   out[1] = prd[1] - 0.5*(prd[2] + prd[3])          #BK
   out[2] = 0.5*(prd[2] + prd[3])                  #BN
   out[3] = 0.25*(prd[2] - prd[3])                 #B±
   out[4] = prd[4]                                #Dab

   out[5] = -(prd[5] + prd[6] + prd[7]) / √3      #T⁰₀(ϵ)
#   out[6] = 0.5*(prd[8] - prd[9])                 #T¹₁(ϵ)
   out[6] = (2*prd[5] - prd[6] - prd[7]) / √6     #T²₀(ϵ)
   out[7] = -prd[8]                               #T²₁(ϵ)
   out[8] = (prd[6] - prd[7])*0.5                 #T²₂(ϵ)
   
   out[ 9] = prd[9]                              #T²₀(χ)
   out[10] = -√(2.0/3.0)*prd[10]                  #T²₁(χ)
   out[11] = prd[11] / √(6.0)                     #T²₂(χ)
   #the factor of 0.5 on the x terms is from 2lₓ = l_+ + l_-
   #out[13] = prd[13]*csl                          #F
   #out[14] = -2.0*prd[13]*prd[14]*csl             #ρzF
   #out[15] = -prd[13]*prd[15]*csl                 #ρxF
   #out[16] = prd[16]*0.5*csl                      #Vn/2
   #out[17] = prd[17]                              #ηz
   #out[18] = 0.5*prd[18]                          #ηx
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
function irrepswap(irrep::String,sdict)
   if irrep=="Il"
      sdict["A"] = 1
      sdict["B"] = 3
      sdict["C"] = 2
   elseif irrep=="IIr"
      sdict["A"] = 2
      sdict["B"] = 3
      sdict["C"] = 1
   elseif irrep=="IIl"
      sdict["A"] = 2
      sdict["B"] = 1
      sdict["C"] = 3
   elseif irrep=="IIIr"
      sdict["A"] = 3
      sdict["B"] = 1
      sdict["C"] = 2
   elseif irrep=="IIIl"
      sdict["A"] = 3
      sdict["B"] = 2
      sdict["C"] = 1
   else #Ir is default
      sdict["A"] = 1
      sdict["B"] = 2
      sdict["C"] = 3
   end
   return sdict
end

tornamind(str::String)::Int = occursin("_",str) ? parse(Int,split(str,"_")[2]) : 1
function torstgsetter(stg::Vector{Int},ctrl)::Vector{Int}
   if ctrl["stages"] > 1
      for i in eachindex(stg)
         stg[i] = stg[i]>0 ? 2 : stg[i]
      end#for
   end#if
end

function secordinp(molnam::String,ctrl)
#The new second order reader will maintain the hard coded forms of:
# H_rot, H_sr, H_qua
#All the torsional parts will be moved to the ℋ vector
#initialize a f v3 rhoz rhox array for each rotor and fill in with read params
#then adjustvalues of rhos by corresponding Fs
#then convert to Op structure and start ℋ
   blk = blockfind(molnam, "%2NDORDER")
   len = blk[2] - blk[1] - 1
   #println(len)
#      ℋ = Op(1.0,tp=[2;0;;]) + Op(0.5) - Op(1.0,tp=[1+iseven(ctrl["NFOLD"][1]);0;;])
#      stgs = [1;1;-1]
   ℋ = Op(0.0)
   stg = [0]
   secns = secordinit_lim()
   secns = irrepswap(ctrl["Irrep"],secns)
   file = readdlm(pwd()*"/"*molnam*".inp",';', skipstart=blk[1])#,comments=true,comment_char='#')
   #println(file)
   tpz = zeros(Int,2,length(ctrl["NFOLD"]))
   val = zeros(Float64,12)
   err = zeros(Float64,12)
   for i in 1:len
      nam = string(strip(file[i,1]))
      #ind = secns[nam]
      if nam ∈ keys(secns)
         ind = secns[nam] #get(secns, nam, 19)
         val[ind] = file[i,2]
         err[ind] = file[i,3]
      elseif occursin("F",nam)
         n = tornamind(nam)
         tmp = zeros(Int,2,length(ctrl["NFOLD"]))
         tmp[1,n] = 2
         ℋ += Op(csl*file[i,2],tp=tmp)
         err = vcat(err,file[i,3])
         stg = vcat(stg,1)
      elseif occursin("V",nam)
         n = tornamind(nam)
         tmp = zeros(Int,2,length(ctrl["NFOLD"]))
         tmp[2,n] = 1+iseven(ctrl["NFOLD"][n])
         ℋ += Op(0.5*csl*file[i,2],tp=zero(tmp)) + Op(-1.0,tp=tmp)
         err = vcat(err,[file[i,3];0])
         stg = vcat(stg,[1;-1])
      elseif occursin("ρ",nam)
         #okay so users will be restriced to an F ρz ρx ordering
         n = tornamind(nam)
         tmp = zeros(Int,2,length(ctrl["NFOLD"]))
         tmp[1,n] = 1
         if occursin("x",nam)
            ℋ += Op(-2*csl*file[i,2]*file[i-2,2],tp=tmp,rp=[1],rf=[nx])
            val[2] += csl*file[i-2,2]*file[i,2]^2 #B -> Beff
         else#this is for z!!!!
            ℋ += Op(-2*csl*file[i,2]*file[i-1,2],d=1,tp=tmp)
            val[1] += csl*file[i-1,2]*file[i,2]^2 #A -> Aeff
            @show file[i-1,2]*file[i,2]
         end
         err = vcat(err,file[i,3])
         stg = vcat(stg,2)
      elseif occursin("η",nam)
         n = tornamind(nam)
         tmp = zeros(Int,2,length(ctrl["NFOLD"]))
         tmp[1,n] = 1
         if occursin("x",nam)
            ℋ += Op(-2*csl*file[i,2]*file[i-2,2],tp=tmp,rp=[1],rf=[sx])
         else
            ℋ += Op(-2*csl*file[i,2]*file[i-2,2],tp=tmp,rp=[1],rf=[sz])
         end
         err = vcat(err,file[i,3])
         stg = vcat(stg,2)
      else
         @warn "Oops! $nam isn't implemented at 2nd order"
      end#else
   end#for
   pop!(val)
   @show val
   @show ℋ
   popat!(err,12)
   popfirst!(ℋ)
   popfirst!(stg)
   val = sod2prep_lim(val)
   #val, err = epszxcheck!(val,err)
   #@show val
   return val, err, ℋ, stg
end

function unitdict()::Dict{String,Float64}
out = Dict{String,Float64}("MHz"=>1.,"cm-1"=>29979.2458,"kHz"=>1e-3,"Hz"=>1e-6,
   "mHz"=>1e-9,"GHz"=>1e3,"THz"=>1e6,"arb"=>1.,
   "eV"=>241_798_840.7662022,"Hart"=>6_579_681_360.732768)
end

function opparse(v::Float64,s::String,list::Dict{String,Op})::Op
   s = split(s,' ')
   out = Op(v)
   for i in eachindex(s)
      part = split(s[i],"^")
      if length(part)==1
         out *= list[s[i]]
      elseif length(part)==2
         out *= list[part[1]]^parse(Int,part[2])
      else
         println("don't raise operators to multiple powers")
      end
   end
   return out
end

function valset(nams,vals,unts,stg)::Vector{Float64}
#this doesn't do what I want it to do. I want a way to comment out a line so
# it isn't used in the current calculation but it is carried over into the 
# new input file. The flaw with this design is the wasted allocations of 
# calculating the matrix elements and then multiplying by 0
#I leave this as an exercise for Sophie
   conv = unitdict()
   for i in eachindex(vals)
      #vals[i] *= nam[i][1]≠'#'
      if stg[i] > 0
         vals[i] *= conv[unts[i]] * (nams[i][1]≠'#')
      #stage-based variable setting will happen in the tsrcalc function 
      #this is so the derivatives know how to handle the added operators
      end#if
   end#for
   return vals
end#
function stgvalset(H,stg)
   for i in eachindex(stg)
      if stg[i] < 0
         H[i].v *= H[i+stg[i]].v
      end
   end
   return H
end

function opreader(molnam,ctrl,ℋ,stgs,errs)
   molnam = "test_input"
   blk = blockfind(molnam, "%OPS")
   blk[1] += 1
   len = blk[2] - blk[1] + 1
   ctrl["opbk"] = blk

   file = try
      readdlm(pwd()*"/"*molnam*".inp",';', skipstart=blk[1],comments=true,
         comment_char='#')
   catch
      zeros(0,0)
   end

   len = size(file,1)

   if size(file,2)!=7
   @warn "You are using a deprecated input format! Check the github for the new structure"
   end

if len != 0 #if there are added parameters
   nams = strip.(file[:,1])
   vals = Float64.(file[:,2])
   unts = string.(strip.(file[:,3]))
   errs = vcat(errs,file[:,4])
   nstgs = Int.(file[:,5])
   vibstring = strip.(file[:,6])
   Opers = string.(strip.(file[:,7]))
   vals = valset(nams,vals,unts,nstgs)

   opdict = Opsdict()
   ℋ += opparse(vals[1],Opers[1],opdict)
   for i in 2:len
      #this is technically concatination but I don't know how to initialize an
      #array of Ops of specified length efficiently
      ℋ += opparse(vals[i],Opers[i],opdict)
   end#for
   #vals, unts, ops are collected in ℋ and vibstr isn't used right now
end#if
   stgs = vcat(stgs,nstgs)
   #ℋ, μs = stgextract(ℋ,stg,0)
   return ℋ,stgs,errs
end#function 

function stgextract(H,stg,n)
   yprm = [stg .== n]
   nprm = [stg .≠ n]
   sub = H[yperm]
   H = H[nperm]
   Hstg = stg[nprm]
   sstg = stg[yprm]
   return H, sub
end

#=
ind = get(baseops,split[1],22)


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
=#
