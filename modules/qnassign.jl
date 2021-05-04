"""
This program contains the qunatum number assignment algorithms

This is a difficult task because each matrix is built from a basis set
    of |J N K m σ> and produces (vtcalc+1)*(4*J+2) eigenvalues with the
    |J N τ vt σ> of these QNs, only J and σ are assigned for free
"""

module qn
using LinearAlgebra
using LinearAlgebra.BLAS
#Utility functions
function closest_index(a, x)
    ibest = 1
    dxbest = abs(a[ibest]-x)
    for i in eachindex(a)
        dx = abs(a[i]-x)
        if dx < dxbest
            dxbest = dx
            ibest = i
        end
    end
    return ibest
end
function Indexer(J,N,τ)
    oldjs = collect(0.5:(J-1))
    index = sum(4 .* oldjs .+ 2)
    index += (N-J+1/2)(2N-1) + N+τ+1
    return index
end
function Δlist(J,S)
    max = J+S
    min = abs(J-S)
    return collect(min:max)
end
function KaSpan(N)
	τs = collect(Int64, -N:N)
	τs = τs .+ N
	kas = floor.(0.5*(τs+ones(size(τs))))
	return kas
end
function KcSpan(N)
	τs = collect(Int64, -N:N)
	τs = τs .+ N
	kc = 2*N .- τs .+ ones(size(τs))
	kc = floor.(0.5*kc)
end
function BigKaSpan(N,Nm)
	N = Int64(N)
	Nm = Int64(Nm)
	if isodd(N)&&isodd(Nm)
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kas = floor.(0.5*(τs+ones(size(τs))))
	elseif iseven(N)&&isodd(Nm)
		Nm += 1
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kas = floor.(0.5*(τs+ones(size(τs))))
		kas = kas[2:end-1]
	elseif isodd(N)&&iseven(Nm)
		Nm += 1
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kas = floor.(0.5*(τs+ones(size(τs))))
		kas = kas[2:end-1]
	elseif iseven(N)&&iseven(Nm)
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kas = floor.(0.5*(τs+ones(size(τs))))
	else
		println("huh, this can't be right in the quantum number assignment")
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kas = floor.(0.5*(τs+ones(size(τs))))
	end
	return kas
end
function BigKcSpan(N,Nm)
	N = Int64(N)
	Nm = Int64(Nm)
	if isodd(N)&&isodd(Nm)
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kc = 2*Nm .- τs .+ ones(size(τs))
		kc = floor.(0.5*kc)
	elseif iseven(N)&&isodd(Nm)
		Nm += 1
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kc = 2*Nm .- τs .+ ones(size(τs))
		kc = floor.(0.5*kc)
		kas = kas[2:end-1]
	elseif isodd(N)&&iseven(Nm)
		Nm += 1
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kc = 2*Nm .- τs .+ ones(size(τs))
		kc = floor.(0.5*kc)
		kas = kas[2:end-1]
	elseif iseven(N)&&iseven(Nm)
		τs = collect(Int64, -N:N)
		τs = τs .+ Nm
		kc = 2*Nm .- τs .+ ones(size(τs))
		kc = floor.(0.5*kc)
	else
		println("huh, this can't be right in the quantum number assignment")
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kas = floor.(0.5*(τs+ones(size(τs))))
	end
	return kas
end

#Assignment Routines
function assign(J::Float64,sigma::Float64,vals,vecs,tvecs)
    Nu = Int64(J+0.5)
    Nl = Int64(J-0.5)
    Nueigs = zeros(2*Nu+1,vtcalc+1)
    Nleigs = zeros(2*Nl+1,vtcalc+1)
    for i in 1:length(vals)
        energy = vals[i]
        vect = vecs[:,i]
        #we are using the highest value in the eigenvector to understand where in
        #the original matrix each eigenvalue came from to figure out what it's QNs must be
        index = iamax(vecs[:,i])
        #this if statement determines if N = J±1/2
        #we should try to assign vt by projecting the K spanning torsional wavefunction
        #    on to a subset of the full wavefunction
        vtcheck = zeros(Float64,vtcalc+1)
        if index<=((vtcalc+1)*(2*J))
            N = Nl
            vtp = convert(Int64,floor(index/(2*N+1.01)))+1
            Nleigs[findfirst(isequal(0.0),Nleigs[:,vtp]),vtp] = energy
        else
            N = Nu
            #Then we use it to determine vt
            index = index-(vtcalc+1)*(2*J)
            vtp = convert(Int64,floor(index/(2*N+1.01)))+1
            Nueigs[findfirst(isequal(0.0),Nueigs[:,vtp]),vtp] = energy
        end
    end
    #sort by energy to "assign" τ
    #This will need to be changed but like not right now
    vmaxind = vtmax+1
    Nueigs = sort(Nueigs[:,1:vmaxind],dims=1)
    Nleigs = sort(Nleigs[:,1:vmaxind],dims=1)
    Jeigs = vcat(Nleigs[:,1:vmaxind],Nueigs[:,1:vmaxind])
    return Jeigs
end

function sortbased(N,vtc,vtm,vals,vecs)
    #I wanted to call this primitive but that's a type name
    Nu = Int64(N)
    Nl = Nu-1
    J = Nu - 0.5
    Nueigs = zeros((2*Nu+1)*(vtc+1))
    Nleigs = zeros((2*Nl+1)*(vtc+1))
    for i in 1:length(vals)
        energy = vals[i]
        vect = vecs[:,i]
        index = iamax(vecs[:,i])
        vtcheck = zeros(Float64,vtc+1)
        if index<=((vtc+1)*(2*J))
            N = Nl
            Nleigs[findfirst(isequal(0.0),Nleigs)] = energy
        else
            N = Nu
            index = index-(2*J)
            Nueigs[findfirst(isequal(0.0),Nueigs)] = energy
        end
    end
    Nueigs = reshape(sort(Nueigs),2*Nu+1,vtc+1)#,dims=1)
    Nleigs = reshape(sort(Nleigs),2*Nl+1,vtc+1)
    vtmaxind = vtm+1
    Jeigs = vcat(Nleigs[:,1:vtmaxind],Nueigs[:,1:vtmaxind])
    return Jeigs
end

function sortbsdv2(N,vtc,vtm,vals,vecs)
    #This new version of the primitive uses subsets of overlap integrals
    #    as opposed to just the largest magnitue element
    Nu = Int64(N)
    Nl = Nu-1
    J = Nu - 0.5
    tJ = 2*Nu-1
    Nueigs = zeros((2*Nu+1)*(vtc+1))
    Nleigs = zeros((2*Nl+1)*(vtc+1))
    for i in 1:length(vals)
        energy = vals[i]
        vect = vecs[:,i]
        #index = iamax(vecs[:,i]) here is the change
        ovNl = transpose(vect[1:(vtc+1)*(tJ)])*vect[1:(vtc+1)*(tJ)]
        ovNu = transpose(vect[(vtc+1)*(tJ)+1:end])*vect[(vtc+1)*(tJ)+1:end]
        #vtcheck = zeros(Float64,vtc+1)
        if ovNl > ovNu
            N = Nl
            Nleigs[findfirst(isequal(0.0),Nleigs)] = energy
        else
            N = Nu
            #index = index-(2*J)
            Nueigs[findfirst(isequal(0.0),Nueigs)] = energy
        end
    end
    Nueigs = reshape(sort(Nueigs),2*Nu+1,vtc+1)#,dims=1)
    Nleigs = reshape(sort(Nleigs),2*Nl+1,vtc+1)
    vtmaxind = vtm+1
    Jeigs = vcat(Nleigs[:,1:vtmaxind],Nueigs[:,1:vtmaxind])
    return Jeigs
end

function suboverlap(Nm,N,vtc,vtm,vals,vecs,σ)
    J = N - 0.5
	#println("J=$J")
	#println(vals)
    Nu = Int64(N)
    Nl = N-1
    Nld = Int64(2*Nl+1)
    Nud = Int64(2*Nu+1)
	tJ = 2*Nu-1
	vtmaxind = vtm+1
    Nueigs = zeros(Float64, (Nud),(vtc+1))
    Nleigs = zeros(Float64, (Nld),(vtc+1))
    Nuvecs = zeros(Float64, (vtc+1)*(Nld+Nud), Nud, (vtc+1))
    Nlvecs = zeros(Float64, (vtc+1)*(Nld+Nud), Nld, (vtc+1))
    for i in 1:length(vals)
        energy = vals[i]
        vect = vecs[:,i]
        offset = (vtc+1)*Nld
        ovNls = zeros(Float64, vtc+1)
        ovNus = zeros(Float64, vtc+1)
        for v in 0:vtc
            ovNls[v+1] = transpose(vect[v*(Nld)+1:(v+1)*(Nld)])*vect[v*(Nld)+1:(v+1)*(Nld)]
			piece = transpose(vect[offset+v*(Nud)+1:offset+(v+1)*(Nud)])
            ovNus[v+1] = piece*vect[offset+v*(Nud)+1:offset+(v+1)*(Nud)]
#ovNus[v+1] = transpose(vect[offset+v*(Nud)+1:offset+(v+1)*(Nud)])
#*vect[offset+v*(Nud)+1:offset+(v+1)*(Nud)]
        end
        lmain = iamax(ovNls)
        umain = iamax(ovNus)
        ucm = ovNus[umain]
        lcm = ovNls[lmain]
#        println("um=$umain,lm=$lmain,ucm=$ucm,lcm=$lcm,energy=$energy")
        #Assignments seem good up to the above println
        if lcm > ucm
            p = findfirst(isequal(0.0),Nleigs[:,lmain])
			#println("p=$p, lmain=$lmain")
            Nleigs[p] = energy
            Nlvecs[:,p,lmain] = vect
        else
            p = findfirst(isequal(0.0),Nueigs[:,umain])
			#println("p=$p, umain=$umain")
            Nueigs[p] = energy
            Nuvecs[:,p,umain] = vect
        end
    end
	#Jeigs is energies sortd by increasing tau and each column is a differnt vt
    #Jvecs each column is a vector for a given tau and going across dim 2 incrases tau, dim 3 is vt
    #Jeigs[l,m] corresponds to Jvecs[:,l,m]
	#QNs is a 2D array with QNs[l,:,m] as the quantum numbers for Jeigs[l,m]
	#2J N Ka Kc v
    Jeigs = vcat(Nleigs[:,1:vtmaxind],Nueigs[:,1:vtmaxind])
	Nlvecs = Nlvecs[:,:,1:vtmaxind]
	Nuvecs = Nuvecs[:,:,1:vtmaxind]
    Jvecs = hcat(Nlvecs,Nuvecs)
#	println("size before padding")
#	println(size(Jvecs))
	#okay so now we need to vertically split the array along the N changes
	#then we will pad with zeros
	vecsnl = Jvecs[1:Nld,:,:]
	vecsnu = Jvecs[Nld+1:end,:,:]
	lpad = Int64((Nm-Nl)*(vtc+1))
	lzpad = zeros(Float64, lpad,Nld+Nud,vtmaxind)
	vecsnl = vcat(lzpad,vecsnl,lzpad)
	upad = Int64((Nm-Nu)*(vtc+1))
	uzpad = zeros(Float64, upad,Nld+Nud,vtmaxind)
	vecsnu = vcat(uzpad,vecsnu,uzpad)
	Jvecs = vcat(vecsnl,vecsnu)
#	println("size after padding")
#	println(size(Jvecs))
	QNs = vcat(fill(Int64(Nl),Nld),fill(Int64(Nu),Nud))
	QNs = hcat(QNs,vcat(KaSpan(Nl),KaSpan(Nu)),vcat(KcSpan(Nl),KcSpan(Nu)),fill(Int64(σ),Nld+Nud))
	QNs = hcat(fill(tJ, Nld+Nud),QNs, zeros(Int64, Nld+Nud))
	for v in 1:vtm
		temp = QNs[:,:,1]
		temp[:,end] .= v
	end
    return Jeigs, Jvecs, QNs
end

function suboverlapv2(J,S,vtc,vtm,vals,vecs)
#=
Okay the goal is to construct a sortperm array so that we can rearrange the
    eigenvalues and the eigenvectors according to the quantum numbers
Set up a range of N values from triangle conditions
Make an array of 0 Ints with Ncount columns that are 2Nmax+1 long
For each eigenvalue, do an overlap integral on each of the N ranges
Assign an N value based on the largest overlap integral
    This will be done by writing the eigval # to the frist 0 in the proper N column
K_a and K_c will be assigned by sorting definitions
The sorting array will be made by taking the first 2N+1 elements of each column
    and making them into a singular array by filter(!iszero,A)
Lastly the eigenvectors will be sorted to match
=#
    Ns = Δ(J,S)
    Ndegns = 2 .* Ns .+ 1
    Neigs = zeros(Float64, Ndegns[end], (vtc+1), length(Ns))
    perms = zeros(Float64, Ndegns[end], (vtc+1), length(Ns))
    for i in 1:length(vals)
        energy = vals[i]
        vect = vecs[:,i]
        offset = zeros(Int64,length(Ns))
        offset[2:end] = (vtc+1).*Ns[1:end-1]
        ovrlp = zeros(Float64,vtc+1,length(Ns))
        for v in 0:vtc
            for n in 1:length(Ns)
				piece = transpose(vect[offset[n]+v*Ndegns[n]+1:offset[n]+(v+1)*Ndegns[n]])
                ovrlp[v+1,n] = abs(piece*vect[offset[n]+v*Ndegns[n]+1:offset[n]+(v+1)*Ndegns[n]])
            end
        end
        #largest intergral is at argmax(ovrlp)
        println(ovrlp)
        Vind = argmax(ovrlp)[1]
        Nind = argmax(ovrlp)[2]
        Neigs[findfirst(isequal(0.0),Neigs[:,Nind,Vind])] = energy
        perms[findfirst(isequal(0.0),perms[:,Nind,Vind])] = i
    end
    vtmaxind = vtm+1
    Neigs = filter(!iszero,Neigs)
    Nvecs = vecs[filter(!iszero,perms)]
    return Neigs, Nvecs
end


#=
function local(J::Float64,sigma::Float64,vals,vecs,vt)
    Nu = Int64(J+0.5)
    Nl = Int64(J-0.5)
    Nueigs = zeros(2*Nu+1)
    Nleigs = zeros(2*Nl+1)
    vtp = vt+1
    for i in 1:length(vals)
        energy = vals[i]
        vect = vecs[:,i]
        index = iamax(vecs[:,i])
        vtcheck = zeros(Float64,vtcalc+1)
        if index<=(2*J)
            N = Nl
            Nleigs[findfirst(isequal(0.0),Nleigs)] = energy
        else
            N = Nu
            index = index-(2*J)
            Nueigs[findfirst(isequal(0.0),Nueigs)] = energy
        end
    end
    Nueigs = sort(Nueigs)#,dims=1)
    Nleigs = sort(Nleigs)#,dims=1)
    Jeigs = vcat(Nleigs,Nueigs)
    return Jeigs
end

function localglobal(N,rprms,sprms,sigma::Float64,gvals,gvecs,tvals,tvecs)
    #so remember that this N is defined as J+1/2
    J = N-0.5
    locals = Array{Float64}(undef,convert(Int64,4*J+2),0)
    #N = convert(Int64,J+0.5)
    for vt in 0:vtmax
        if J==0.5
            locmat = Hlocalf0(rprms,sprms,vt,sigma,tvals,tvecs)
        else
            locmat = Hlocal(N,rprms,sprms,vt,sigma,tvals,tvecs)
        end
        locval,locvec = eigen!(locmat)
        locval = qnlocal(J,sigma,locval,locvec,vt)
#        println("J=$J")
#        println(locval)
#        println(gvals)
        locals = hcat(locals,locval)
    end
    #now that we have all of our local values, we need to compare them to the global ones
    yikes = Array{Int64}(undef,0,1)
    outvals = zeros(Float64,size(locals))
    vecsize = convert(Int64,(vtcalc+1)*(4*J+2))
    outvecs = zeros(Float64,vecsize,vecsize,vtmax+1)
    #vecord = zeros(Int64,4*J+2,vtmax+1) #This will be used to resort the vectors
    for vtp in 1:(vtmax+1)
        for i in 1:convert(Int64,4*J+2)
            #This iterates down the list of the J block
            ind = closest_index(gvals,locals[i,vtp])
            outvals[i,vtp] = gvals[ind]
            outvecs[:,i,vtp] = gvecs[:,ind]
        end
    end
    if -V3^2 in locals
        println("Should have seen this coming... J=$J, σ=$sigma")
    else
        nothing
    end
    #I need to add an eigenvector return as wells
    return outvals[:,1:vtmax+1]#,outvecs[:,:,1:vtmax+1]
end
=#
function simfact(J::Float64,N,gvals,gvecs,ovals,truncvecs)
    if N < J
        #first check No=N, Jo=J-1
        oldstart = Indexer(J-1,N,-N)
        oldend = Indexer(J-1,N,N)
        prevecs = truncvecs[:,oldstart:oldend]
        SFs = transpose(prevecs)*gvecs
        bestSF = BLAS.iamax(SFs)
    else
    end
end


end

#now i need a function to diagonalize and assign quantum numbers
#find index of highest value in eigenvector, assign QNs as function of index
#I'm going to just use bins and then definitionally assign Ks
#the bins will be temorary 2D arrays, one for each N
#then each column will share an vt value
#finally each column will be sorted to assign the K
#=
function qnassign(J::Float64,sigma::Float64,vals,vecs,tvecs)
    Nu = Int64(J+0.5)
    Nl = Int64(J-0.5)
    Nueigs = zeros(2*Nu+1,vtcalc+1)
    Nleigs = zeros(2*Nl+1,vtcalc+1)
    for i in 1:length(vals)
        energy = vals[i]
        vect = vecs[:,i]
        #we are using the highest value in the eigenvector to understand where in
        #the original matrix each eigenvalue came from to figure out what it's QNs must be
        index = iamax(vecs[:,i])
        #this if statement determines if N = J±1/2
        #we should try to assign vt by projecting the K spanning torsional wavefunction
        #    on to a subset of the full wavefunction
        vtcheck = zeros(Float64,vtcalc+1)
        if index<=((vtcalc+1)*(2*J))
            N = Nl
            vtp = convert(Int64,floor(index/(2*N+1.01)))+1
            Nleigs[findfirst(isequal(0.0),Nleigs[:,vtp]),vtp] = energy
        else
            N = Nu
            #Then we use it to determine vt
            index = index-(vtcalc+1)*(2*J)
            vtp = convert(Int64,floor(index/(2*N+1.01)))+1
            Nueigs[findfirst(isequal(0.0),Nueigs[:,vtp]),vtp] = energy
        end
    end
    #sort by energy to "assign" τ
    #This will need to be changed but like not right now
    vmaxind = vtmax+1
    Nueigs = sort(Nueigs[:,1:vmaxind],dims=1)
    Nleigs = sort(Nleigs[:,1:vmaxind],dims=1)
    Jeigs = vcat(Nleigs[:,1:vmaxind],Nueigs[:,1:vmaxind])
    return Jeigs
end

function qnlocal(J::Float64,sigma::Float64,vals,vecs,vt)
    Nu = Int64(J+0.5)
    Nl = Int64(J-0.5)
    Nueigs = zeros(2*Nu+1)
    Nleigs = zeros(2*Nl+1)
    vtp = vt+1
    for i in 1:length(vals)
        energy = vals[i]
        vect = vecs[:,i]
        index = iamax(vecs[:,i])
        vtcheck = zeros(Float64,vtcalc+1)
        if index<=(2*J)
            N = Nl
            Nleigs[findfirst(isequal(0.0),Nleigs)] = energy
        else
            N = Nu
            index = index-(2*J)
            Nueigs[findfirst(isequal(0.0),Nueigs)] = energy
        end
    end
    Nueigs = sort(Nueigs)#,dims=1)
    Nleigs = sort(Nleigs)#,dims=1)
    Jeigs = vcat(Nleigs,Nueigs)
    return Jeigs
end

function qnlocalglobal(N,rprms,sprms,sigma::Float64,gvals,gvecs,tvals,tvecs)
    #so remember that this N is defined as J+1/2
    J = N-0.5
    locals = Array{Float64}(undef,convert(Int64,4*J+2),0)
    #N = convert(Int64,J+0.5)
    for vt in 0:vtmax
        if J==0.5
            locmat = Hlocalf0(rprms,sprms,vt,sigma,tvals,tvecs)
        else
            locmat = Hlocal(N,rprms,sprms,vt,sigma,tvals,tvecs)
        end
        locval,locvec = eigen!(locmat)
        locval = qnlocal(J,sigma,locval,locvec,vt)
        println("J=$J")
        println(locval)
#        println(gvals)
        locals = hcat(locals,locval)
    end
    #now that we have all of our local values, we need to compare them to the global ones
    yikes = Array{Int64}(undef,0,1)
    outvals = zeros(Float64,size(locals))
    vecsize = convert(Int64,(vtcalc+1)*(4*J+2))
    outvecs = zeros(Float64,vecsize,vecsize,vtmax+1)
    #vecord = zeros(Int64,4*J+2,vtmax+1) #This will be used to resort the vectors
    for vtp in 1:(vtmax+1)
        for i in 1:convert(Int64,4*J+2)
            #This iterates down the list of the J block
            #ind = findnearest(gvals,locals[i,vtp])[1]
            ind = closest_index(gvals,locals[i,vtp])
            outvals[i,vtp] = gvals[ind]
            outvecs[:,i,vtp] = gvecs[:,ind]
        end
    end
    #println(outvals)
    if -V3^2 in locals
        println("Should have seen this coming... J=$J, σ=$sigma")
    else
        nothing
    end
    #I need to add an eigenvector return as wells
    return outvals[:,1:vtmax+1]#,outvecs[:,:,1:vtmax+1]
end

function qnsimfact(J::Float64,N,gvals,gvecs,ovals,truncvecs)
    if N < J
        #first check No=N, Jo=J-1
        oldstart = Indexer(J-1,N,-N)
        oldend = Indexer(J-1,N,N)
        prevecs = truncvecs[:,oldstart:oldend]
        SFs = transpose(prevecs)*gvecs
        bestSF = BLAS.iamax(SFs)
    else
    end
end
=#
