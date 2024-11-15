
function unitdict()::Dict{String,Float64}
out = Dict{String,Float64}("MHz"=>1.,"cm-1"=>29979.2458,"kHz"=>1e-3,"Hz"=>1e-6,
   "mHz"=>1e-9,"GHz"=>1e3,"THz"=>1e6)
end


function Opparse(s::String,list::Dict{String,Op};v=1.0)::Op
   s = split(s,' ')
   out = Op(v)
   for i in eachindex(s)
      part = split(s,"^")
      if length(part)==1
         out *= Op(1.0,list[s[1]])
      elseif length(part)==2
         out *= Op(1.0,list[s[1]])^parse(Int,s[2])
      else
         println("don't raise Operators to multiple powers")
      end
   end
   return out
end
function valset(vals,unts,stg)::Vector{Float64}
   conv = unitdict()
   for i in eachindex(vals)
      if stg[i] > 0
      vals[i] *= conv[unts[i]]
      end#if
   end#for
   return vals
end#

function valset_withcomments(nams,vals,unts,stg)::Vector{Float64}
#this doesn't do what I want it to do. I want a way to comment out a line so
# it isn't used in the current calculation but it is carried over into the 
# new input file. I leave this as an exercise for SOphie
   conv = unitdict()
   for i in eachindex(vals)
      if stg[i] > 0
      vals[i] *= conv[unts[i]] * (1 - (nam[i][1]=='#'))
      end#if
   end#for
   return vals
end#

function Opreader(molnam)
   molnam = "test_input"
   blk = blockfind(molnam, "%OPS")
   blk[1] += 1
   len = blk[2] - blk[1] + 1

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

if len != zero(len) #if there are added parameters
   nams = strip.(file[:,1])
   vals = file[:,2]
   unts = strip.(file[:,3])
   errs = file[:,4]
   stgs = file[:,5]
   vibstring = strip.(file[:,6])
   Opers = strip.(file[:,7])
   vals = valset(vals,unts,stgs)
   
   ℋ = Opparse(vals[1],Opers[1])
   for i in 2:len
      #this is technically concatination but I don't know how to initialize an
      #array of Ops of specified length efficiently
      ℋ += Opparse(vals[i],Opers[i])
   end#for
   return ℋ,
end#function 


#=
ind = get(baseOps,split[1],22)


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
