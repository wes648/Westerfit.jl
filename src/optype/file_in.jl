
function unitdict()::Dict{String,Float64}
out = Dict{String,Float64}("MHz"=>1.,"cm-1"=>29979.2458,"kHz"=>1e-3,"Hz"=>1e-6,
   "mHz"=>1e-9,"GHz"=>1e3,"THz"=>1e6,"nm"=>2.99792458e11,"arb"=>1.,
   "eV"=>241_798_840.7662022,"Hart"=>6_579_681_360.732768)
end

function opparse(s::String,list::Dict{String,Op};v=1.0)::Op
   s = split(s,' ')
   out = Op(v)
   for i in eachindex(s)
      part = split(s,"^")
      if length(part)==1
         out *= Op(1.0,list[s[1]])
      elseif length(part)==2
         out *= Op(1.0,list[s[1]])^parse(Int,s[2])
      else
         println("don't raise operators to multiple powers")
      end
   end
   return out
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
#The new second order reader will maintain the hard coded forms of:
# H_rot, H_sr, H_qua
#All the torsional parts will be moved to the ℋ vector
#initialize a f v3 rhoz rhox array for each rotor and fill in with read params
#then adjustvalues of rhos by corresponding Fs
#then convert to Op structure and start ℋ
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
      vals[i] *= conv[unts[i]] * nam[i][1]≠'#'
      end#if
   end#for
   return vals
end#
function opreader(molnam,ctrl)
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
   vals = file[:,2]
   unts = strip.(file[:,3])
   errs = file[:,4]
   stgs = file[:,5]
   vibstring = strip.(file[:,6])
   Opers = strip.(file[:,7])
   vals = valset(nams,vals,unts,stgs)
   
   ℋ = opparse(vals[1],Opers[1])
   for i in 2:len
      #this is technically concatination but I don't know how to initialize an
      #array of Ops of specified length efficiently
      ℋ += opparse(vals[i],Opers[i])
   end#for
   #vals, unts, ops are collected in ℋ and vibstr isn't used right now
   return ℋ,stgs,errs
else
   return nothing,nothing,nothing
end#if
end#function 


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
