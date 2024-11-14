molnam = "test_input"
 blk = blockfind(molnam, "%OPS")
 blk[1] += 1
 len = blk[2] - blk[1] + 1

   file = try
      readdlm(pwd()*"/"*molnam*".inp",';', skipstart=blk[1],comments=true,comment_char='#')
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
   unts = file[:,3]
   errs = file[:,4]
   stgs = file[:,5]
   vibstring = strip.(file[:,6])
   opers = strip.(file[:,7])

   for i in 1:len
   		entry = opers[i]
   		split = split(entry," ")
   	    for i in 1:length(split)
   	    	if occursin("^",split[i]) == true
   	    		splitspecific = split[i]
   	    		loc = findfirst("^",splitspecific)[1]
   	    		justop = splitspecific[1:loc-1]
   	    		num = parse(Int64,splitspecific[loc+1:end])
   				ind = get(baseops,justop,22)
   				
   			else
   				ind = get(baseops,split[i],22)
   				val[ind] = vals[i]
   				err[ind] = errs[i]
   			end



baseops = Dict("Nz" => 1, "N2" => 2, "Np" => 3, "Nm" => 4, "Npm" => 5, "Nx" => 6,
 "sNy"=> 7, "NS" => 8, "S2" => 9, "Sz" => 10, "Sp" => 11, "Sm" => 12, "Spm" => 13,
 "Sx"=> 14, "sSy" => 15, "Pα" => 16, "cosα" => 17, "Pβ" => 18, "cosβ" => 19,
  "Pγ" => 20, "cosγ" => 21)
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

