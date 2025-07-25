
@kwdef struct Controls
   NFOLD::Vector{Int} = [0]
   S::Float64 = 0.0
   TK::Float64 = 8.0
   mcalc::Int = 8
   vtmax::Int = 0
   Jmax::Float64 = 0.
   apology::Bool = true
   νmin::Float64 = 0.0
   νmax::Float64 = 40.0
   INTthres::Float64 = 0.0001
   λlm0::Float64 = 0.001
   RUNmode::String = "ESF"
   turducken::Int = 1
   maxiter::Int = 60
   overwrite::Bool = true 
   assign::String = "ram36"
   REJECT::Float64 = 10.0
   Irrep::String = "Ir"
   goal::Float64 = 1.0
   mmax::Int = 6
   stages::Int = 1 
   ctbk::Vector{Int} = zeros(Int,2)
   sobk::Vector{Int} = zeros(Int,2)
   inbk::Vector{Int} = zeros(Int,2)
   opbk::Vector{Int} = zeros(Int,2)
end

function ctrlinit()
   ctrl = Dict("NFOLD" => 0, "S" => 0., "TK" => 8.0, "mcalc" => 8, "vtmax" => 0,
      "Jmax" => 0, "apology" => true, "νmin"=>0.0, "νmax"=>40., "INTthres"=>0.00001, 
      "λlm0"=>0.0001, "RUNmode"=>"ESF", "turducken"=>1, "maxiter"=>60, "overwrite"=>true,
      "assign"=>"ram36", "REJECT"=>1.0e+1, "Irrep"=>"Ir", "goal"=>1.0, "mmax"=>6, "stages"=>1,
      "ctbk"=>[0;0],"sobk"=>[0;0],"inbk"=>[0;0],"opbk"=>[0;0])
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
   @show blknam
   @show out
   return out
end
function blockfind_all(molnam::String)#,ctrl)
   count = 0
   inblock = false
   blck = "init"
   blknams = ["%CNTRLS"; "%2NDORDER"; "%INTS"; "%OPS"]
   ctrlnms = ["ctbk"; "sobk"; "inbk"; "opbk"]
   ctrl = ctrlinit()
   io = open("$molnam.inp")
   for i in eachline(io)
      count += 1
      check = contains.(i,blknams)
      if true ∈ check
         blck = ctrlnms[findfirst(check)]
         ctrl[blck][1] = count
         inblock = true
      elseif isempty( filter(x -> !isspace(x), i) )&& inblock==true
         ctrl[blck][2] = count
         #break
      elseif eof(io)
         ctrl[blck][2] = count
         break
      end
   end
   close(io)
   ctrl["opbk"][1] +=1
   #@show ctrl.ctbk 
   #@show ctrl.sobk 
   #@show ctrl.inbk 
   #@show ctrl.opbk 
   return ctrl
end

function ctrlinp(molnam::String,dctrl)
#   ctrl = ctrlinit()
   blk = dctrl["ctbk"]  #blockfind(molnam, "%CNTRLS")
   len = blk[2] - blk[1] - 1
#   ctrl.ctbk  = blk
   file = readdlm(pwd()*"/"*molnam*".inp",'=', skipstart=blk[1],comments=true,comment_char='#')
   for i in 1:len
      nam = strip(file[i,1])
      val = file[i,2]
      if (nam == "NFOLD")&&typeof(val)==Int
         val = [val]
      end
      if typeof(val)==String||typeof(val)==SubString{String}
         val = String(strip(val))
      end
      dctrl[nam] = val
   end
   ctrl = Controls(; (Symbol.(keys(dctrl)) .=> values(dctrl))... )
   #trlS  = Float64(ctrl.S )
   if ctrl.apology  == true
      println("Sorry about the name...")
   end
   if iseven(2*ctrl.Jmax )&&isodd(2*ctrl.S )
      ctrl.Jmax  += 0.5
      println("Jmax must be half integer for half interger S. Adding 1/2")
   elseif isodd(2*ctrl.Jmax )&&iseven(2*ctrl.S )
      ctrl.Jmax  += 0.5
      println("Jmax must be integer for interger S. Adding 1/2")
   end
   if ctrl.NFOLD ==0
      ctrl.mcalc  = 0
      ctrl.vtmax  = 0
   elseif length(ctrl.NFOLD ) > 1
      println("\nYou are running westerfit in ntop mode")
      println("By using this mode, you agree to not complain about the runtime of the program\n")
   end
   #ctrl.assign  = strip(ctrl.assign )
   @show ctrl.assign
   #ctrl.Irrep  = String(strip(ctrl.Irrep ))
   #println(ctrl)
   return ctrl
end

function secordinit_lim()::Dict{String,Int}
   prd = Dict("A" => 1, "B" => 2, "C" => 3, "Dab" => 4,
              "Z" => 1, "X" => 2, "Y" => 3, "Dxz" => 4, 
      "ϵzz" => 5, "ϵxx" => 6, "ϵyy" => 7, "ϵzx" => 8, "ϵxz" => 8,
      "Czz" => 5, "Cxx" => 6, "Cyy" => 7, "Czx" => 8, "Cxz" => 8,
      "χzz"=> 9, "χxz"=> 10, "χxmy"=> 11, "χxx-χyy"=>11,
      "α"=>9, "δ"=>10, "β"=>11)#,
      #"F" => 13, "ρz" => 14, "ρ" => 14, "ρx" =>15, "Vn" => 16, "V3" =>16,
      #"ηz" => 17, "η" => 17, "ηx" => 18)
   return prd
end
function secnam_init()::Vector{String}
   return ["AZ"; "BX"; "CY"; "Dab"; "ϵzz"; "ϵxx"; "ϵyy"; "ϵzx"; "χzz"; "χxx-χyy"; "χxz"]
end
function sod2prep_lim(prd::Array{Float64},nf)::Array{Float64}
   out = zeros(11)
   ##tempa = prd[1] + csl*prd[13]*prd[14]^2         #Aeff = A + Fρz²
   ##tempb = prd[2] + csl*prd[13]*prd[15]^2         #Beff = B + Fρx²
   out[1] = prd[1] - 0.5*(prd[2] + prd[3])        #BK
   out[2] = 0.5*(prd[2] + prd[3])                 #BN
   out[3] = 0.25*(prd[2] - prd[3])                #B±
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
   for i ∈ 1:length(nf)
      out[ 9+4i] = prd[13]*csl                          #F
      out[10+4i] = -2.0*prd[13]*prd[14]*csl             #ρzF
      out[11+4i] = -prd[13]*prd[15]*csl                 #ρxF
      out[12+4i] = prd[16]*0.5*csl                      #Vn/2
   end
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
      sdictA  = 1
      sdictB  = 3
      sdictC  = 2
   elseif irrep=="IIr"
      sdictA  = 2
      sdictB  = 3
      sdictC  = 1
   elseif irrep=="IIl"
      sdictA  = 2
      sdictB  = 1
      sdictC  = 3
   elseif irrep=="IIIr"
      sdictA  = 3
      sdictB  = 1
      sdictC  = 2
   elseif irrep=="IIIl"
      sdictA  = 3
      sdictB  = 2
      sdictC  = 1
   else #Ir is default
      sdictA  = 1
      sdictB  = 2
      sdictC  = 3
   end
   return sdict
end

tornamind(str::String)::Int = occursin("_",str) ? parse(Int,split(str,"_")[2]) : 1
function torstgsetter(stg::Vector{Int},ctrl)::Vector{Int}
   if ctrl.stages  > 1
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
   blk = ctrl.sobk  #blockfind(molnam, "%2NDORDER")
   len = blk[2] - blk[1] - 1
   #println(len)
#      ℋ = Op(1.0,tp=[2;0;;]) + Op(0.5) - Op(1.0,tp=[1+iseven(ctrl.NFOLD [1]);0;;])
#      stgs = [1;1;-1]
   ℋ = Vector{Op}(undef,0)
   nams = secnam_init()
   stg = zeros(Int,0)
   secns = secordinit_lim()
   secns = irrepswap(ctrl.Irrep ,secns)
   file = readdlm(pwd()*"/"*molnam*".inp",';', skipstart=blk[1])#,comments=true,comment_char='#')
   #println(file)
   tpz = zeros(Int,2,length(ctrl.NFOLD ))
   val = zeros(Float64,11+4*length(ctrl.nf))
   errs = zeros(Float64,length(val))
   unts = fill("MHz",11)
   for i in 1:len
      nam = string(strip(file[i,1]))
      nvl = file[i,2]
      #ind = secns[nam]
      if nam ∈ keys(secns)
         ind = secns[nam] #get(secns, nam, 19)
         val[ind] = nvl
         errs[ind] = file[i,3]
      elseif occursin("F",nam) && !iszero(nvl)
         n = tornamind(nam)
         tmp = zeros(Int,2,length(ctrl.NFOLD ))
         tmp[1,n] = 2
         push!(val,csl*nvl)
         push!(ℋ, Op(nam,tp=tmp))
         push!(errs,file[i,3])
         push!(stg,1)
         push!(unts,"MHz")
      elseif occursin("V",nam) && !iszero(nvl)
         n = tornamind(nam)
         tmp = zeros(Int,2,length(ctrl.NFOLD ))
         tmp[2,n] = 1+iseven(ctrl.NFOLD[n])
         push!(ℋ, Op(nam,rf=[OpFunc(E,1)],tp=zero(tmp)), Op(nam,tp=tmp))
         push!(val, 0.5*csl*nvl, -1.0)
         push!(errs,file[i,3],0)
         push!(stg,1,-1)
         push!(unts,"cm-1","cm-1")
      elseif occursin("ρ",nam) && !iszero(nvl)
         #okay so users will be restriced to an F ρz ρx ordering
         n = tornamind(nam)
         tmp = zeros(Int,2,length(ctrl.NFOLD ))
         tmp[1,n] = 1
         if occursin("x",nam)
            push!(ℋ, Op(nam,tp=tmp,rf=[OpFunc(nx,1)] ))
            push!(val,-2*csl*nvl*file[i-2,2])
            val[2] += csl*file[i-2,2]*nvl^2 #B -> Beff
         else#this is for z!!!!
            push!(ℋ, Op(nam,d=1,tp=tmp))
            push!(val,-2*csl*nvl*file[i-1,2])
            val[1] += csl*file[i-1,2]*nvl^2 #A -> Aeff
            #@show file[i-1,2]*file[i,2]
         end
         push!(errs,file[i,3])
         push!(stg,0)
         push!(unts,"unitless")
      elseif occursin("η",nam) && !iszero(nvl)
         n = tornamind(nam)
         tmp = zeros(Int,2,length(ctrl.NFOLD ))
         tmp[1,n] = 1
         if occursin("x",nam)
            push!(ℋ, Op(nam,tp=tmp,rf=[OpFunc(sx,1)] ))
         else
            push!(ℋ, Op(nam,tp=tmp,rf=[OpFunc(sz,1)] ))
         end
         push!(val,nvl)
         push!(errs,file[i,3])
         push!(stg,0)
         push!(unts,"MHz")
      elseif iszero(nvl)
      else
         @warn "Oops! $nam isn't implemented at 2nd order"
      end#else
   end#for
   val[1:11] = sod2prep_lim(val[1:11])
   #val, err = epszxcheck!(val,err)
   #@show val
   return val, errs, ℋ, stg, unts
end

function unitdict()::Dict{String,Float64}
out = Dict{String,Float64}("MHz"=>1.,"cm-1"=>29979.2458,"kHz"=>1e-3,"Hz"=>1e-6,
   "mHz"=>1e-9,"GHz"=>1e3,"THz"=>1e6,"arb"=>1.,"z"=>0.0,
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

function opstrproc!(T::Type,abcd::Vector{Int},tops::Array{Int,2},fvec::Vector{OpFunc},
                    str::String,pow::Int,dict::Dict{String,Function}
                    )::Tuple{Vector{Int},Array{Int,2},Vector{OpFunc}}
   rcheck = findall(isequal(str),["N2";"NS";"S2";"Nz"] )
   tcheck = findall(isequal(str),["Pα" "Pβ" "Pγ"; "cosα" "cosβ" "cosγ"] )
   if isempty(rcheck) && isempty(tcheck)
      push!( fvec, OpFunc(T,dict[str],pow) )
   elseif !isempty(rcheck) && isempty(tcheck)
      abcd[rcheck] .+= pow
   elseif isempty(rcheck) && !isempty(tcheck)
      tops[tcheck] .+= pow
   else
      @warn "This operator isn't defined. Contact westerfit support for help"
      println("Seriously Wes needs more people to talk to about this")
   end
   return abcd,tops,fvec
end
function opparse2(nam::String,s::String,lnf::Int,list::Dict{String,Function})::Op
   s = string.(split(s,' '))
   fvec = Vector{OpFunc}(undef,0)
   partype = Float64
   typecheck = 1
   abcd = zeros(Int,4)
   tops = zeros(Int,2,lnf)
   for i in eachindex(s)
      part = string.(split(s[i],"^"))
      #if part[1] ∈ imlist
      #   typecheck *= -1
      #   partype = ComplexF64
      #end
      if length(part)==1
         abcd,tops,fvec = opstrproc!(partype,abcd,tops,fvec,s[i],1,list)
      elseif length(part)==2
         abcd,tops,fvec = opstrproc!(partype,abcd,tops,fvec,part[1],parse(Int,part[2]),list)
      else
         println("don't raise operators to multiple powers")
      end
      partype = Float64
   end
   out = Op(nam,fvec,abcd,tops)
   return out
end

function valset(nams,vals,unts,stg)::Vector{Float64}
#this doesn't do what I want it to do. I want a way to comment out a line so
# it isn't used in the current calculation but it is carried over into the 
# new input file. The flaw with this design is the wasted allocations of 
# calculating the matrix elements and then multiplying by 0
#I leave this as an exercise for Sophie
   conv = unitdict()
   for i in eachindex(nams)
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

function opreader(molnam,ctrl,vals,errs,ℋ,stgs)
   #molnam = "test_input"
   blk = ctrl.opbk  #blockfind(molnam, "%OPS")
   #blk[1] += 1
   len = blk[2] - blk[1] + 1
   #ctrl.opbk  = blk
   lnf = length(ctrl.NFOLD )
   file = try
      readdlm(pwd()*"/"*molnam*".inp",';', skipstart=blk[1],comments=true,
         comment_char='#')
   catch
      zeros(0,0)
   end

   len = size(file,1)

   if size(file,2)≠7
   @warn "You are using a deprecated input format! Check the github for the new structure"
   end

if len ≠ 0 #if there are added parameters
   nams = string.(strip.(file[:,1]))
   append!(vals,Float64.(file[:,2]))
   unts = string.(strip.(file[:,3]))
   append!(errs,file[:,4])
   nstgs = Int.(file[:,5])
   vibstring = strip.(file[:,6])
   Opers = string.(strip.(file[:,7]))
   vals = valset(nams,vals,unts,nstgs)

   opdict = Opsdict()
   push!(ℋ,opparse2(nams[1],Opers[1],lnf,opdict))
   for i in 2:len
      #this is technically concatination but I don't know how to initialize an
      #array of Ops of specified length efficiently
      #ℋ += opparse(vals[i],Opers[i],opdict)
      push!(ℋ,opparse2(nams[i],Opers[i],lnf,opdict))
   end#for
   #vals, unts, ops are collected in ℋ and vibstr isn't used right now
end#if
   append!(stgs,nstgs)
   #ℋ, μs = stgextract(ℋ,stg,0)
   return ℋ,stgs,errs,unts
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

function μparse(val::Float64,op::String,dict)
   op = split.(split(op), "^")
   out = μOb(val,
      μFuncR(dict[op[1][1]],parse(Int,op[1][2]),0),
      μFuncT(dict[op[2][1]],parse(Int,op[2][2]),0))
   return out
end

function intreader(molnam,ctrl)
   blk = ctrl.inbk  #blockfind(molnam, "%OPS")
   @show blk
   blk[1] += 1
   len = blk[2] - blk[1] - 1
   ctrl.inbk  = blk
   lnf = length(ctrl.NFOLD )
   file = try
      readdlm(pwd()*"/"*molnam*".inp",';', skipstart=blk[1],comments=true,
         comment_char='#')
   catch
      zeros(0,0)
   end
   #len = size(file,1)
   #if size(file,2)≠7
   #@warn "You are using a deprecated input format! Check the github for the new structure"
   #end
   @show len
   @show size(file)
   file = file[1:len,:]
   @show size(file)
if len ≠ 0 #if there are added parameters
   nams = string.(strip.(file[:,1]))
   @show file[:,2]
   vals = Float64.(file[:,2])
   #units is 3
   #vibs is 4
   ops = string.(strip.(file[:,5]))
   μs = Vector{μOb}(undef,0)
   intdict = Intdict()
   push!(μs, μparse(vals[1],ops[1], intdict) )
   for i ∈ 2:len
      push!(μs, μparse(vals[i],ops[i], intdict) )
   end
else
   @warn "Yikes! Need some dipoles for a simulator"
end#if
   return μs
end#function


function lineprep(lns,nf,s,vtm)
   #converts the input file into a more code friendly format
   #           1  2  3   4   5  6  7  8   9  10  11   12
   #input1 = [ju nu kau kcu mu jl nl kal kcl ml freq unc]
   #           1  2  3   4   5   6  7  8  9   10  11 12  13   14
   #input2 = [ju nu kau kcu vtu σu jl nl kal kcl vtl σl freq unc]
   #           1  2   3   4   5   6
   #output = [ju scu indu jl scl indl]
if length(lns[1,:]) == 12
   nf = iseven(nf) ? Int(0.5*nf) : nf
   qunus = lns[:,1:10]
   freqs = lns[:,11]
   uncs = lns[:,12]
   inds = zeros(Int,size(lns,1),6)
   #J_u
   inds[:,1] = Int.(2 .* qunus[:,1])
   #ind_u
   inds[:,3] = qn2ind.(nf,vtm,qunus[:,5],qunus[:,1],s,qunus[:,2],qunus[:,3],qunus[:,4])
   #J_l
   inds[:,4] = Int.(2 .* qunus[:,6])
   #ind_l
   inds[:,6] = qn2ind.(nf,vtm,qunus[:,10],qunus[:,6],s,qunus[:,7],qunus[:,8],qunus[:,9])
   if !iszero(nf) #derive sc from m
      #sc_u
      inds[:,2] = Int.(mod.(qunus[:, 5],nf)) .+ 1
      #sc_l
      inds[:,5] = Int.(mod.(qunus[:,10],nf)) .+ 1
   else
      inds[:,2] .= 1
      inds[:,5] .= 1
   end
elseif length(lns[1,:]) == 14 #still working on
   nf = iseven(nf) ? Int(0.5*nf) : nf
   qunus = lns[:,1:12]
   freqs = lns[:,13]
   uncs = lns[:,14]
   qunus[:,  5] = @. ceil(Int,0.5*lns[:, 5])*powneg1(isodd(lns[:, 5]))*nf
   qunus[:, 11] = @. ceil(Int,0.5*lns[:,11])*powneg1(isodd(lns[:,11]))*nf

   inds = zeros(Int,size(lns,1),6)
   inds[:,1] = Int.(2 .* qunus[:,1])
   inds[:,3] = qn2ind.(nf,vtm,qunus[:,5],qunus[:,1],s,qunus[:,2],qunus[:,3],qunus[:,4])
   inds[:,4] = Int.(2 .* qunus[:,7])
   inds[:,6] = qn2ind.(nf,vtm,qunus[:,11],qunus[:,7],s,qunus[:,8],qunus[:,9],qunus[:,10])

   inds[:,2] = Int.(mod.(qunus[:, 6],nf)) .+ 1
   inds[:,5] = Int.(mod.(qunus[:,11],nf)) .+ 1

else #length(lns[1,:]) != 12
    println("Wrong number of columns in .lne file")
end
   return inds, freqs, uncs
end

function qn2ind(nf,vtm,m,j,s,n,ka,kc)
   if (nf≠zero(nf))
   #ka = abs(ka)
   ind = sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)*(vtm+1)
   #ind += (vtm+floor(Int,m/nf))*(2*j+1) #<--- This function is wrong for vt≠0
   ind *= Int(2s+1)
   ind += (ceil(Int,0.5*vtm) + floor(Int,m/nf))*(2*j+1)
   ind += sum(2 .* collect((j-s):(n-1)) .+ 1) + n + kakc2k(n,ka,kc) + 1
   ind = convert(Int,ind)
   return ind
   else
   return qn2ind(j,s,n,ka,kc)
   end
end
function qn2ind(j,s,n,ka,kc)
   jp = (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)
   np = sum(2 .* collect((j-s):(n-1)) .+ 1)
   kp = n + kakc2k(n,ka,kc) + 1
   ind = jp + np + kp
   ind = convert(Int,ind)
   return ind
end
function kakc2k(n,ka,kc)
   ka = abs(ka)
   if iseven(n+ka+kc)
      return ka*powneg1(n+ka)
   else isodd(n+ka+kc)
      return -ka*powneg1(n+ka)
   end
#   return ka*powneg1(n+ka+kc)
end
