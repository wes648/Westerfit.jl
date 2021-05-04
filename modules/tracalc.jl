
module tracalc
using DelimitedFiles
using Printf
using LinearAlgebra
include("./WIGXJPF.jl")
using .WIGXJPF

"""
This calculates the transitions from the energy levels. It is called TraCalc
Basis set: Ψ=|J N τ vt σ⟩

Potential improvements:
	speed
	remove use of hardcoded index based selection rules

Index = Σ_{j=1/2}^{J-1}(4j+2) + (N-J+1/2)(2N-1) + N+τ+1
Transitions will be done by explicit expressions with the 6 possibilities
ΔJ= 0,ΔN=+1,ΔK=-1  a
    index=2J
    (N+2-K)(N+1-K)/(4N+4)

ΔJ= 0,ΔN=+1,ΔK= 0  b
    index=2J+1
    (N+1+K)(N+1-K)/(N+1)

ΔJ=+1,ΔN= 0,ΔK=-1  a
    index=2N+2
    (N+1-K)(N+K)(2N+1)/4N(N+1)

ΔJ=+1,ΔN= 0,ΔK= 0  b
    index=2N+3
    (2N+1)K^2/(N^2+N)

ΔJ=+1,ΔN=+1,ΔK=-1  a
    index=2J+2N+5
    (N+2-K)(N+1-K)/(4N+4)

ΔJ=+1,ΔN=+1,ΔK= 0  b
    index=2J+2N+6
    (N+1+K)(N+1-K)/(N+1)

We will only calculate the ΔJ+ΔN≥0 for the speed. This is justified by the program
    not being intended for vibrational spectra
We can replace the Wigner symbols with direct expressions. Wolfram has the 6-j I need
The rest of the 3-js should be easy to find
I also need to make a list of the excluded terms for the a-types as a function of
    Nmax. Nmax=2 excludes 2,5
    arr = collect(1:Nmax-1)
    arr = collect(Iterators.flatten(zip(arr,arr)))
    exclude = zeros(length(arr))
    exclude[1] = 2
    for i=1:length(arr)
        exclude[i+1] = exclude[i]+2*arr[i]+1
    end
"""

#TOP='PROLATE'
μa = 1.0
μb = 1.0
#molnam="BELGI"

function T1mu(q)
	if q==0
		μ=μa
	elseif q==1
		μ=μb
	elseif q==-1
		μ=μb
	else
		println("This dipole can't be right")
		μ=μa
	end
	return μ
end


function Δlist(J,S)
    max = J+S
    min = abs(J-S)
    return collect(min:max)
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
		kas = kas[1:end-2]
	elseif isodd(N)&&iseven(Nm)
		Nm += 1
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kas = floor.(0.5*(τs+ones(size(τs))))
		kas = kas[1:end-2]
	elseif iseven(N)&&iseven(Nm)
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kas = floor.(0.5*(τs+ones(size(τs))))
	else
		println("huh, this can't be right...")
		τs = collect(Int64, -Nm:Nm)
		τs = τs .+ Nm
		kas = floor.(0.5*(τs+ones(size(τs))))
	end
	return kas
end


function TraCalchCs(molnam,Nmax,energies,sigma,vt)
    #much of the following array creation should be moved into a separate function as TraCalc
    #This is because it will be the same for each vt and σ pair.
	#Doing this first could open up parallel options
    arr = collect(Int64,1:Nmax-1)
    arr = collect(Int64,Iterators.flatten(zip(arr,arr)))
    #this excluded list is for the levels that won't have a-type transitions
    #this isn't needed at all
    #exclude = zeros(Int64,length(arr))
    #exclude[1] = 2
    #for i=1:(length(arr)-1)
    #    exclude[i+1] = exclude[i]+2*arr[i]+1
    #end
    #this list has the J value for the level of identical index
    posjs = collect(Float64,0.5:Nmax-0.5)
#    println(posjs)
    jinds = fill(posjs[1],Int64(4*posjs[1]+2))
    for i = 2:length(posjs)
        jtemp = fill(posjs[i],Int64(4*posjs[i]+2))
        jinds = vcat(jinds,jtemp)
    end
    #jtemp = fill(posjs[end],Int64(2*posjs[end]))
    #jinds = vcat(jinds,jtemp)
#    println(jinds)
    #this list has the N value for the level of indentical index
    ninds = zeros(Int64,length(energies))
    #this list has the K value for the level of indentical index
    kinds = zeros(Int64,length(energies))
    shift=1
    for i=1:length(posjs)
        inds = shift+1
        #indm = shift+2+2*i
        indf = 4*i+2+shift
        ninds[inds:indf] = fill(i,indf-inds+1)
        newks = collect(Int64,-i:i)
        newks = vcat(newks,newks)
        kinds[inds:indf] = newks
        shift += 4*i+2
    end
    ninds = ninds[1:length(jinds)]
    kinds = kinds[1:length(jinds)]
	#freqs is actually frequencies[1] and intensities[2] and type[3]
    freqs = Array{Float64}(undef,0,3)
	#states contains the indices of the states for a given frequency
    states = Array{Int64}(undef,0,2)
    leng = size(energies)[1]
    #there has to be a better way to do the following but alas
    sigind = sigma + 1
    vtind = vt + 1
    for i=1:leng
        N = ninds[i]
        K = kinds[i]
        J = jinds[i]
        #ΔJ= 0,ΔN=+1,ΔK= 0  b
        indexu = i+2*J+1
        indexu = convert(Int64,indexu)
        if (indexu≤leng)&&(J-N>0)
            freq = energies[indexu,vtind,sigind]-energies[i,vtind,sigind]
            intent = (μb^2)*(N+1+K)*(N+1-K)/(N+1)
            type = 1.0
            if freq > 0
#                freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [indexu i]
                states = vcat(states,state)
            else
                freq = abs(freq)
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [i indexu]
                states = vcat(states,state)
            end
        else
            nothing
        end
        #ΔJ=+1,ΔN= 0,ΔK= 0  b
        indexu = i+2*N+1
        indexu = convert(Int64,indexu)
        if (indexu≤leng)&&(J-N<0)
            freq = energies[indexu,vtind,sigind]-energies[i,vtind,sigind]
            intent = (μb^2)*(2*N+1)*K^2/(N^2+N)
            type = 2.0
            if freq > 0
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [indexu i]
                states = vcat(states,state)
            else
                freq = abs(freq)
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [i indexu]
                states = vcat(states,state)
            end
        else
            nothing
        end
        #ΔJ=+1,ΔN=+1,ΔK= 0  b
        indexu = i+2*J+2*N+4
        indexu = convert(Int64,indexu)
        if indexu≤leng
            freq = energies[indexu,vtind,sigind]-energies[i,vtind,sigind]
            intent = (μb^2)*(N+1+K)*(N+1-K)/(N+1)
            type = 3.0
            if freq > 0
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [indexu i]
                states = vcat(states,state)
            else
                freq = abs(freq)
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [i indexu]
                states = vcat(states,state)
            end
        else
            nothing
        end
        #calcas
        #ΔJ= 0,ΔN=+1,ΔK=-1  a
        indexu = i+2*J
        indexu = convert(Int64,indexu)
        if (indexu≤leng)&&(J-N>0)
            freq = energies[indexu,vtind,sigind]-energies[i,vtind,sigind]
            intent = (μa^2)*(N+2-K)*(N+1-K)/(4*N+4)
            type = 4.0
            if freq > 0
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [indexu i]
                states = vcat(states,state)
            else
                freq = abs(freq)
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [i indexu]
                states = vcat(states,state)
            end
        else
            nothing
        end
        #ΔJ=+1,ΔN= 0,ΔK=-1  a
        indexu = i+2*N+2
        indexu = convert(Int64,indexu)
        if (indexu≤leng)&&(J-N<0)
            freq = energies[indexu,vtind,sigind]-energies[i,vtind,sigind]
            intent = (μa^2)*(N+1-K)*(N+K)*(2*N+1)/(4*N*(N+1))
            type = 5.0
            if freq > 0
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [indexu i]
                states = vcat(states,state)
            else
                freq = abs(freq)
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [i indexu]
                states = vcat(states,state)
            end
        else
            nothing
        end
        #ΔJ=+1,ΔN=+1,ΔK=-1  a
        indexu = i + 2*J+2*N+3
        indexu = convert(Int64,indexu)
        if indexu≤leng
            freq = energies[indexu,vtind,sigind]-energies[i,vtind,sigind]
            intent = (μa^2)*(N+2-K)*(N+1-K)/(4*N+4)
            type = 6.0
            if freq > 0
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [indexu i]
                states = vcat(states,state)
            else
                freq = abs(freq)
                #freqline = [freq intent]
                freqline = [freq intent type]
                freqs = vcat(freqs,freqline)
                state = [i indexu]
                states = vcat(states,state)
            end
        else
            nothing
        end
    end
    #At this point there are two vital arrays: freqs and states
    neworder = sortperm(freqs[:,1])
    freqs = freqs[neworder,:]
    states = states[neworder,:]
    if sigma==0
        ju = jinds[states[:,1]]
        ju = ju+fill(0.5,size(ju))
        jl = jinds[states[:,2]]
        jl = jl+fill(0.5,size(jl))
        nu = ninds[states[:,1]]
        nl = ninds[states[:,2]]

		tauinds = kinds + ninds
		ka = floor.(0.5*(tauinds+ones(size(kinds))))
        kc = 2*ninds-tauinds+ones(size(kinds))
        kc = floor.(0.5*kc)

		ku = kinds[states[:,1]]
        ku += nu
        kau = ku+ones(size(ku))
        kau = floor.(0.5*(kau))
        kcu = 2*nu+ones(size(ku))-ku
        kcu = floor.(0.5*kcu)
        kl = kinds[states[:,2]]
        kl += nl
        kal = kl+ones(size(kl))
        kal = floor.(0.5*(kal))
        kcl = 2*nl+ones(size(kl))-kl
        kcl = floor.(0.5*kcl)
        sigs = fill(sigma, length(ju))
        open("AStates_$molnam.txt", "w") do io
            writedlm(io, [nu kau kcu ju  nl kal kcl jl sigs sigs freqs])
#            writedlm(io, [ju nu ku sigs jl nl kl sigs freqs])
        end
		egys = energies[:,vtind,sigind]#./29979.2458
		open("AStates_$molnam.egy","w") do io
			writedlm(io, [egys ninds ka kc jinds], ',')
		end
    else
        ju = jinds[states[:,1]]
        jl = jinds[states[:,2]]
        nu = ninds[states[:,1]]
        nl = ninds[states[:,2]]
        ku = kinds[states[:,1]]
		tauinds = kinds + ninds
		ka = floor.(0.5*(tauinds+ones(size(kinds))))
        kc = 2*ninds-tauinds+ones(size(kinds))
        kc = floor.(0.5*kc)
		kau = floor.(0.5*(ku+ones(size(ku))))
        kcu = 2*nu+ones(size(ku))-ku
        kcu = floor.(0.5*kcu)
        kl = kinds[states[:,2]]
        kal = floor.(0.5*(kl+ones(size(kl))))
        kcl = 2*nu+ones(size(kl))-kl
        kcl = floor.(0.5*kcl)
        sigs = fill(sigma, length(ju))
        open("EStates_$molnam.txt", "w") do io
            writedlm(io, [ju nu kau kcu sigs jl nl kal kcl sigs freqs])
#            writedlm(io, [ju nu ku sigs jl nl kl sigs freqs])
        end
		egys = energies[:,vtind,sigind]#./29979.2458
		open("EStates_$molnam.egy","w") do io
			writedlm(io, [egys ninds ka kc jinds])
		end
    end
    return freqs, states
end

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

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

#=function RotInt(Nu,uvec,Nl,lvec)
	kau = KaSpan(Nu)
	kal = KaSpan(Nl)
	out = 0.0
	for i in 1:length(kau)
		part = uvec[i]*(-1)^(Nl-kau[[i]-1])
		for j in 1:length(kal)
			part *= lvec[j]
			for q in -1:1
				out += part*wig3j(Nl,1,Nu,Kl,q,-Ku)*T1mu(q)
			end#closes q loop
		end#closes kl loop
	end#closes ku loop
	out *= (2.0*Nu+1.0)*sqrt(2.0*Nl+1.0)
end=#

function RotInt(Nu,uvec,Nl,lvec)
	kau = collect(-Nmax:Nmax)
	out = 0.0
	for i in 1:length(kau)
		part = uvec[i]*(-1)^(Nl-kau[[i]-1])
		for q in -1:1
			Kl = Ku - q
			if Nmax>abs(Kl)
				nothing
			else
				j = i - q
				out += part*lvec[j]*wig3j(Nl,1,Nu,Kl,q,-Ku)*T1mu(q)
			end
		end#closes q loop
	end#closes ku loop
	out *= (2.0*Nu+1.0)*sqrt(2.0*Nl+1.0)
end
function SRInt(Ju,Jl,vecu,vecl)
	Nus = Δlist(Ju,S)
	Nls = Δlist(Jl,S)
	vecu = reshape(vecu,(2*Nmax+1,2*S+1))
	vecl = reshape(vecl,(2*Nmax+1,2*S+1))
	out = 0.0
	for i in 1:size(vecu)[2]
		nu = Nus[i]
		for j in 1:size(vecl)[2]
			nl = Nls[j]
			part = RotInt(nu,uvec,nl,lvec)
			part *= sqrt(2*nl+1)*(2*nu+1)*(-1)^(nl+S+Ju+1)*wig6j(nl,Jl,S,Ju,nu,1)
		end#closes nl loop
		out += part
	end#closes nu loop
	out = abs(out)^2
	return out
end

#=function SRInt(Ju,Jl,vecu,vecl)
	Nus = Δlist(Ju,S)
	Nls = Δlist(Jl,S)
	vecu = reshape(vecu,(2*Nmax+1,2*S+1))
	vecl = reshape(vecl,(2*Nmax+1,2*S+1))
	out = 0.0
	for i in 1:size(vecu)[2]
		nu = Nus[i]
		Kaus = BigKaSpan(nu,Nmax)
		for j in 1:size(vecl)[2]
			nl = Nls[j]
			Kals = BigKaSpan(nl,Nmax)
			for k in 1:size(vecu)[1]
				ku = Kaus[k]
				for l in 1:size(vecl)[1]
					kl = Kals[l]
					part = 0.0
					for q in -1:1
						part += wig3j(nl,1,nu,kl,q,-ku)*T1mu(q)
					end#closes q loop
					part *= vecl[l,j]
				end#closes kl loop
				part *= vecl[k,i]*(-1)^(nl-1-ku)
			end#closes ku loop
			part *= sqrt(2*nl+1)*(2*nu+1)*(-1)^(nl+S+Ju+1)*wig6j(nl,Jl,S,Ju,nu,1)
		end#closes nl loop
		out += part
	end#closes nu loop
	out = abs(out)^2
	return out
end=#

function TorInt(Nmax,tvecs,vtu,vtl,ku,kl)
	out = transpose(tvecs[:,vtu+1,Int64(ku)+Nmax+1])*tvecs[:,vtl+1,Int64(kl)+Nmax+1]
	return out
end

function TSRInt(Nmax,S,vtc,Ju,Jl,vtu,vtl,vecu,vecl,tvecs)
	Nus = Δlist(Ju,S)
	Nls = Δlist(Jl,S)
	Sd = Int64(2*S+1)
	Nmd = Int64(2*Nmax+1)*(vtc+1)
	vecu = reshape(vecu,(Nmd,Sd))
	vecl = reshape(vecl,(Nmd,Sd))
	out = 0.0
	part = 0.0
	for i in 1:size(vecu)[2]
		nu = Nus[i]
		#Kaus = repeat(BigKaSpan(nu,Nmax),outer=[vtc+1])
		Kaus = collect(Int64, -nu:nu)
		for j in 1:size(vecl)[2]
			nl = Nls[j]
			#Kals = repeat(BigKaSpan(nl,Nmax),outer=[vtc+1])
			Kals = collect(Int64, -nl:nl)
			for k in 1:size(Kaus)[1]
				ku = Kaus[k]
				for l in 1:size(Kals)[1]
					kl = Kals[l]
					for q in -1:1
						part += wig3j(nl,1,nu,kl,q,-ku)*T1mu(q)*TorInt(Nmax,tvecs,vtu,vtl,ku,kl)
					end#closes q loop
					part *= vecl[l,j]
				end#closes kl loop
				part *= vecl[k,i]*(-1)^(nl-1-ku)
			end#closes ku loop
			part *= sqrt(2*nl+1)*(2*nu+1)*(-1)^(nl+S+Ju+1)*wig6j(nl,Jl,S,Ju,nu,1)
		end#closes nl loop
		out += part
		part = 0.0
	end#closes nu loop
	out = abs(out)^2
	return out
end

function boltzmod(int::Float64, temp::Float64, elower::Float64)::Float64
	out = int*exp(-elower/(temp*2.083661912E+4))
	return out
end

function TraCalcCs(Nmax,S,vtc,energies,vecs,tvecs,EQNs,sigma,vt)
	maxΔJ = 1
	freqmin = 26500
	freqmax = 40000
	intmin = 1E-6
	Tk = 20.0
	freqs = Array{Float64}(undef,0,3) #freq int
	fquns = Array{Int64}(undef,0,12) #2Ju Nu Kau Kcu sigma vt 2Jl Nl Kal Kcl sigma vt
	for i in 1:length(energies)
		for j in i:length(energies)
			if abs(EQNs[i,1,vt+1]-EQNs[j,1,vt+1]) > maxΔJ*2
				break
			else
				freq = energies[j] - energies[i]
				#println("j=$j")
				if freq > 0
					int = TSRInt(Nmax,S,vtc,EQNs[j,1,vt+1]/2,EQNs[i,1,vt+1]/2,vt,
					vt,vecs[:,j],vecs[:,i],tvecs)
					int = boltzmod(int,Tk,energies[i])
					freqline = [freq int energies[i]]
#					println(EQNs[j,:,vt+1])
					qnline = transpose(vcat(EQNs[j,:,vt+1], EQNs[i,:,vt+1]))
#					println(qnline)
#					println(transpose(qnline))
				elseif freq<0
					freq = abs(freq)
					int = TSRInt(Nmax,S,vtc,EQNs[i,1,vt+1]/2,EQNs[j,1,vt+1]/2,vt,
					vt,vecs[:,i],vecs[:,j],tvecs)
					int = boltzmod(int,Tk,energies[j])
					freqline = [freq int energies[j]]
					qnline = transpose(vcat(EQNs[i,:,vt+1], EQNs[j,:,vt+1]))
				else
					nothing
				end#closes if statement
				if (freq<freqmax)&&(freq>freqmin)&&(int > intmin)
					freqs = vcat(freqs,freqline)
					fquns = vcat(fquns,qnline)
				else
					nothing
				end
			end
		end#closes j iterator
	end#closes primary iterator
	fquns = Int64.(fquns)
	return freqs, fquns
end

function TraWriter(molnam,freqs, qunus) #emulates the cat file structure of SPCAT
	c = 29979.2458
	p = sortperm(freqs[:,1])
	freqs = freqs[p,:]
	qunus = qunus[p,:]
	out = fill("0",size(freqs)[1])
	for i in 1:size(freqs)[1]
		#freq
#		part = lpad(@sprintf("%0.4f", freqs[i,1]), 13)
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
		#J N Ka Kc sigma vt is the order in the array
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







end
