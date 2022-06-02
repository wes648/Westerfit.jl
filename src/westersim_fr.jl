"""
Welcome to westersim_fr.jl

This program is a new, more group-theoretically driven version of westersim.
The fr stands for free rotor as this version does not use the belgi 2-stage
    approach and instead does everything in a single stage. The matrix is twice
    Wang transformed to improve eigenvector definition

input & fit parameters:
      1  2  3  4   5   6   7    8    9    10  11
prs =[A; B; C; F; tθ; V3; ϵzz; ϵxx; ϵyy; ϵxzb; η;
        12   13  14  15  16
		ΔN; ΔNK; ΔK; δN; δK;
		17  18  19   20  21  22  23   24   25   26   27   28   29    30   31  32   33
		Fm; V6; V3m; ρm; ρ3; FN; FK; Fbc; Fab; V3N; V3K; V3ab; V3bc; ρN; ρK; ρab; ρbN
		 34   35    36    37    38   39
		ΔsN; ΔsNK; ΔsKN; ΔsK; δsN; δsK;
		40   41   42  43  44  45   46
		ΦJ; ΦJK; ΦKJ; ΦK; ϕJ; ϕJK; ϕK]

internal calculation parameters:
      1  2  3   4    5  6   7   8  9 10 11 12
prs =[A; B; C; Dab; Fr; ρ; V3; ao; a; b; d; η;
		13   14  15  16  17
		ΔN; ΔNK; ΔK; δN; δK;
		18  19   20  21  22  23  24   25   26   27   28   29    30   31  32   33   34
		Fm; V6; V3m; ρm; ρ3; FN; FK; Fbc; Fab; V3N; V3K; V3ab; V3bc; ρN; ρK; ρab; ρbN
		 35   36    37    38   39   40
		ΔsN; ΔsNK; ΔsKN; ΔsK; δsN; δsK;
		41   42   43  44  45   46  47
		ΦJ; ΦJK; ΦKJ; ΦK; ϕJ; ϕJK; ϕK]

"""

using LinearAlgebra
using LinearAlgebra.BLAS
using LinearAlgebra.LAPACK
using DelimitedFiles
using Optim
using SparseArrays
include("/home/wes/files/westerfit/experimental/WIGXJPF.jl")
using .WIGXJPF

apology = true

println("westerfit!")
print("Enter molecule name: \n")
#molnam = readline()
#molnam = "mepho-m"
#include("$molnam.par")
#println(molnam)
#println("Nmax=$Nmax")
if apology==true
	println("Sorry about the name...")
end

global NFOLD=3
global TK = 25.0

################################################################################
############                       westersim                        ############
################################################################################

function Δlist(J,S)
	max = Int64(J+S)
	min = Int64(abs(J-S))
	return collect(min:max)
end
function paramshift(input)
	out = zeros(Float,length(input)+1)
	out[1:3] = input[1:3]
	out[4] = (input[1] - input[2])*input[5] #Dab
	out[5] = input[4]/(input[4]-input[1]) #reduction of F
	out[6] = input[1]/input[4] # ρ
	out[7] = input[7]
	out[8] = -(input[7] + out[9] + out[10])/3.0 #ao
	out[9] = -(2.0*input[7] - out[9] - out[8])/6.0 #a
	out[10] = (input[8] - input[9])*0.5 #b
	out[11] = -input[10] #d
	out[12:end] = input[11:end]
	return out
end
function klist(n)
	Γ1 = sort([collect(-1:-2:-n); collect(0:2:n)]) .+ (n + 1)
	Γ2 = sort([collect(-2:-2:-n); collect(1:2:n)]) .+ (n + 1)
	if iseven(n)
		return Γ1, Γ2
	else
		return Γ2, Γ1
	end
end
function klist(n,m)
	Γ1 = sort([collect(-1:-2:-n); collect(0:2:n)]) .+ (n + 1)
	Γ2 = sort([collect(-2:-2:-n); collect(1:2:n)]) .+ (n + 1)
	if iseven(n)&&(m≥0)||isodd(n)&&(m<0)
		return Γ1, Γ2
	else
		return Γ2, Γ1
	end
end
function mklist(n,mcalc)
	A1, A2 = klist(n,-mcalc)
	nd = convert(Int,2*n+1)
	shift = nd
	for m ∈ (1 - mcalc):mcalc
		a1, a2 = klist(n,m)
		A1 = vcat(A1,a1 .+ shift)
		A2 = vcat(A2,a2 .+ shift)
		shift += nd
	end
	return A1, A2
end
function eye(x)
	spdiagm(ones(x))
end
function θ(j::Float64,n,s::Float64)::Float64
	if s==0.5
		out = (n-j)/(j+0.5)
	#elseif n==0.0
    #    out = 0.0
    else
        out = n*(n+1.0) + s*(s+1.0) - j*(j+1.0)
        out = out/(2.0*n*(n+1.0))
    end
    return out
end
function ϕ(j::Float64,n,s::Float64)::Float64
    if s==0.5
        out = -1.0/(j+0.5)
    else
        out = (n-j+s)*(n+j+s+1)*(s+j-n+1)*(n+j-s)
        out *= 1.0/((2.0*n-1.0)*(2.0*n+1.0))
        out = -sqrt(out)/n
    end
    return out
end
function f(x,y::Float64)::Float64
    out = sqrt(x*(x+1.0) - y*(y+1.0))
    return out
end
function g(x::Float64,y::Float64)::Float64
    out = sqrt((x-y)*(x-y-1.0))
    return out
end
function Δtest(a,b,c)
	return c ⊆ collect(abs(a-b):abs(a+b))
end
#These are the various Wang
function un(m)
	if iseven(Int64(m))
		n = Int64(m/2)
		out = (1/sqrt(2)) .* [-eye(n) rotl90(eye(n)); rotl90(eye(n)) eye(n)]
		return out
	elseif isodd(Int64(m))
		n = Int64((m-1)/2)
		out = (1/sqrt(2)) .* [-eye(n) zeros(n) rotl90(eye(n)); zeros(1,n) sqrt(2) zeros(1,n);
		 rotl90(eye(n)) zeros(n) eye(n)]
		return out
	end
end
function ut(m,n)
   nd = 2*n+1
   out = (1/sqrt(2)) .*  [-1*eye(m) zeros(m) rotl90(eye(m)); zeros(1,m) sqrt(2) zeros(1,m);
    rotl90(eye(m)) zeros(m) eye(m)]
   out = kron(out,eye(nd))
   return out
end
function ut(m,j,s)
	nlist = Δlist(j,s)
	out = zeros(0,0)
	for i in 1:length(nlist)
		out = cat(out,ut(m,nlist[i]);dims=(1,2))
	end
	return out
end
function ur(n,m)
   md = 2*m+1
   out = (1/sqrt(2)) .* [-eye(n) zeros(n) rotl90(eye(n)); zeros(1,n) sqrt(2) zeros(1,n);
    	rotl90(eye(n)) zeros(n) eye(n)]
   out = kron(eye(md),out)
   return out
end
function ur(j,s,m)
	nlist = Δlist(j,s)
	out = zeros(0,0)
	for i in 1:length(nlist)
		out = cat(out,ur(nlist[i],m);dims=(1,2))
	end
	return out
end

function k2ka(n,ka)
	#ka = abs(k)
	return ka
end
function k2kc(n,k)
	ka = abs(k)
	if k < 0
		kc = n - ka + 1
	else
		kc = n - ka
	end
	return kc
end
function qngen(n,m,σ)
	nd = 2*n+1
	md = 2*m+1
	narray = fill(n,nd*md)
	karray = kron(ones(md),collect(-n:n))
	marray = NFOLD .* kron(collect(-m:m),ones(nd)) .+ σ
	σarray = fill(σ,nd*md)
	out = hcat(narray,karray,marray,σarray)
end
function qngen(j,s,m,σ)
	nlist = Δlist(j,s)
	out = zeros(0,4)
	for i in 1:length(nlist)
		out = vcat(out,qngen(nlist[i],m,σ))
	end
	jlist = fill(j,size(out)[1])
	out = hcat(jlist,out)
	return out
end
function qngen2(n,m,σ)
	nd = Int(2*n+1)
	md = Int(2*m+1)
	narray = fill(n,nd*md)
	karray = kron(ones(Int,md),collect(Int,-n:n))
	kaarray = k2ka.(narray,karray)
	kcarray = k2kc.(narray,karray)
	marray = NFOLD .* kron(collect(Int,-m:m),ones(Int,nd)) .+ σ
	σarray = fill(σ,nd*md)
	out = hcat(narray,kaarray,kcarray,marray,σarray)
end
function qngen2(j,s,m,σ)
	nlist = Δlist(j,s)
	out = zeros(0,5)
	for i in 1:length(nlist)
		out = vcat(out,qngen2(nlist[i],m,σ))
	end
	jlist = fill(j,size(out)[1])
	out = hcat(jlist,out)
	return out
end

#ρa = λa*A/F
#ρb = λb*B/F
#r = 1 - (λa^2)*A/F - (λb^2)*B/F
#F = F0/r
function Htorm0(pr,N,K,m,σ)
	out = @. pr[5]*(m+σ)^2  - 2*pr[5]*pr[6]*(m+σ)*K + pr[5]*(pr[6]*K)^2 + pr[7]*0.5 + pr[17]*(m+σ)^2 + 0.5*pr[18]
	out = @. out + pr[20]*K*(m+σ)^3 + pr[23]*(m+σ)^2*K^2 + pr[31]*(m+σ)*K^3
	out = @. out + 2.0*ρ*m*K + pr[27]*K^2 + 2.0*pr[19]*m^2
	out = @. out + pr[22]*N*(N+1.0)*m^2 + pr[30]*m*K*N*(N+1.0) + pr[26]*N*(N+1.0)
end
function Htorm1(pr,N,K,m,σ)
	#<m+1 K|H|m K>
    out = @. -pr[7]*0.25 - pr[21]*K*(m+1.0) - 0.5*pr[19]*(m^2 + (m+1.0)^2)
	out = @. out - 0.5*pr[26]*N*(N+1.0) - 0.5*pr[27]*K^2
end
function Htorm2(pr,N,K,m,σ)
	#<m+2 K|H|m K>
    out = @. -pr[18]*0.25 - 0.0*m*K
end
function Hr0K(pr,N,K)
	out = @. 0.5*(B + pr[3])*N*(N+1.0) + (pr[1] - 0.5*(pr[2] + pr[3]))*K^2
end
function Hr1K(pr,N,K)
	#<N K+1 | H |N K>
	out = @. (pr[4] + pr[28])*sqrt(N*(N+1.0)-K*(K-1.0))*(K-0.5)
end
function Hr2K(pr,N,K)
	out = @. 0.25*(pr[2] - pr[3])  + pr[29]
	out = @. out*sqrt((N*(N+1.0) - K*(K - 1.0))*(N*(N+1.0) - (K - 1.0)*(K - 2.0)))
end
function Hrt1K1m(pr,N,K,m,σ)
	#<N K+1 m+1 | H |N K>
	out = @. -0.5*pr[28]*sqrt(N*(N+1.0)-K*(K-1.0))*(K-0.5) + 0.0*m*K
end
function Hrt2K1m(pr,N,K,m,σ)
	out = @. -0.5*pr[29] + 0.0*m*K
	out = @. out*sqrt((N*(N+1.0) - K*(K - 1.0))*(N*(N+1.0) - (K - 1.0)*(K - 2.0)))
end
function Hrtoffboth(pr,N,m,σ)
	if Float64(N)==0.0
		return [0.0]
	else
		karray = collect(Float64,-N:N)
		nd = length(karray)
		of1diag = Hrt1K1m(pr,N,karray[2:end],m,σ)
		of2diag = Hrt2K1m(pr,N,karray[3:end],m,σ)
		mat = spdiagm(nd,nd,1=>of1diag,2=>of2diag,-1=>of1diag,-2=>of2diag)
#		mat = Symmetric(mat)
		return mat
	end
end
function Hrot(pr,N)
	if Float64(N)==0.0
		return [0.0]
	else
    	karray = collect(Float64,-N:N)
    	nd = length(karray)
    	ondiags = Hr0K(pr,N,karray)
    	of1diag = Hr1K(pr,N,karray[2:end])
    	of2diag = Hr2K(pr,N,karray[3:end])
    	Rotmat = spdiagm(nd,nd,0=>ondiags,1=>of1diag,2=>of2diag,-1=>of1diag,-2=>of2diag)
#    	Rotmat = Symmetric(Rotmat)
    	return Rotmat
	end
end
function Htor(mpr,calc,N,σ)
	tnp = convert(Int,2*N+1)
	if mcalc==0
		return spzeros(tnp,tnp)
	else
	tmp = convert(Int,2*mcalc+1)
	ks = kron(ones(Float64,2*mcalc+1),collect(Float64,-N:N))
	ms = kron(3.0 .* collect(-mcalc:mcalc),ones(Float64,tnp))
	ondiags = Htorm0(pr,N,ks,ms,σ)
	of1diags = Htorm1(pr,N,ks[tnp+1:end],ms[tnp+1:end],σ)
	of2diags = Htorm2(pr,N,ks[tnp+2:end],ms[tnp+2:end],σ)
	out = spdiagm(tnp*tmp,tnp*tmp,0=>ondiags,tnp=>of1diags,tnp+1=>of2diags,-tnp=>of1diags,-tnp-1=>of2diags)
	return out
	end
end
function Htr(pr,N,m,σ)
	nd = 2*N+1
	md = 2*m+1
	out = zeros(Float64,md*nd,md*nd)
	out = kron(eye(md),Hrot(pr,N))
	out += Htor(m,N,σ)
	out += kron(spdiagm(1=>ones(md-1),-1=>ones(md-1)), Hrtoffboth(pr,N,m,σ))
#	return Symmetric(out)
	return out
end
function Hs0K0N(pr,J,N,S,thet,K)
#   <J S N K| H |J S N K>
    out  = @. -0.5*pr[8]*(J*(J+1.0)-N*(N+1)-S*(S+1.0)) + thet*(pr[9]*(3.0*K^2 - N*(N+1)))
end
function Hs1K0N(pr,J,N,thet,K)
#   <J S N K-1| H |J S N K>
    out = @. pr[11]*(K-0.5)*f(N,K-1.0)*thet
end
function Hs2K0N(pr,J,N,thet,K)
    #   <J S N K-2| H |J S N K>
    out = @. thet*0.5*pr[10]*sqrt((N*(N+1)-K*(K-1.0))*(N*(N+1)-(K-1.0)*(K-2.0)))
end
function Hsm2K1N(pr,N::Float64,K::Array{Float64,1},ϕ)::Array{Float64,1}
#   <J S N+1 K-2| H |J S N K>
    #out = @. -0.25*(srp[3]+srp[12]*(K*(N-K)+(K-2)*(N-K+2)))*f(N,K-1)*g(N,-K+1.0)
    out = @. 0.25*pr[10]*f(N+1.0,K-2.0)*g(N+1.0,K-1.0).*ϕ
end
function Hsm1K1N(pr,N::Float64,K::Array{Float64,1},ϕ)::Array{Float64,1}
#   <J S N+1 K-1| H |J S N K>
    out = @. 0.25*pr[11]*(N+2.0*K)*g(N+1.0,K-1.0).*ϕ
end
function Hs0K1N(pr,N::Float64,K::Array{Float64,1},ϕ)::Array{Float64,1}
    #out = @. (1.5*srp[2]*K+0.5*K*(srp[8]*K^2+srp[9]*N^2))*sqrt(N^2 - K^2)
    out = @. 1.5*pr[9]*K*sqrt((N+1.0)^2 - K^2).*ϕ
end
function Hs1K1N(pr,N::Float64,K::Array{Float64,1},ϕ)::Array{Float64,1}
    out = @. 0.25*pr[11]*(N-2.0*K)*g(N+1.0,-K-1.0).*ϕ
end
function Hs2K1N(pr,N::Float64,K::Array{Float64,1},ϕ)::Array{Float64,1}
    #out = @. 0.25*(srp[3]+srp[12]*(K*(N+K)+(K+2)*(N+K+2)))*f(N,K)*g(N,K+1.0)
    out = @. -0.25*pr[10]*f(N+1.0,K+1.0)*g(N+1.0,-K-1.0).*ϕ
end
#Matrix Builder
function Hspi0N(pr,J,S,N)
    thet = θ(J,N,S)
    Nt2 = N*(N+1.0)
    if N==0.0
        #ondiags = Hs0K0N(J,N,S,Nt2,thet,[0.0],srprms)
        Spimat = zeros(Float64,1,1)#diagm(0=>ondiags)
    else
        karray = collect(Float64,-N:N)
        ondiags = Hs0K0N(pr,J,N,S,thet,karray)
        of1diag = Hs1K0N(pr,J,N,thet,karray[2:end])
        of2diag = Hs2K0N(pr,J,N,thet,karray[3:end])
        Spimat = spdiagm(length(karray),length(karray),0=>ondiags,1=>of1diag,2=>of2diag,-1=>of1diag,-2=>of2diag)
    end
#	return Symmetric(Spimat)
	return Spimat
end
function Hspi1N(pr,J::Float64,S::Float64,Nl)
    jϕ = ϕ(J,Nl+1.0,S)
    karray = collect(Float64,-Nl:Nl)
	pm1 = Hsm2K1N(pr,Nl,karray[2:end],jϕ)
	p0 = Hsm1K1N(pr,Nl,karray,jϕ)
	p1 = Hs0K1N(pr,Nl,karray,jϕ)
	p2 = Hs1K1N(pr,Nl,karray,jϕ)
	p3 = Hs2K1N(pr,Nl,karray[1:end-1],jϕ)
    Spimat = spdiagm(2*Int(Nl)+1,2*Int(Nl)+3, -1=>pm1, 0=>p0, 1=>p1, 2=>p2, 3=>p3)
    return Spimat
end
function srprep(J,S)
	ns = Δlist(J,S)
	nd = 2 .* Int.(ns) .+ 1
	out = ones(Int, length(ns),2)
	out[1,2] = nd[1]
	for i in 2:length(ns)
		out[i,1] = out[i-1,2] + 1
		out[i,2] = out[i,1] + nd[i] - 1
	end
	jd = Int((2*S+1)*(2*J+1))
	return ns, nd, out, jd
end
function srprep(J,S,md)
	ns = Δlist(J,S)
	nd = 2 .* Int.(ns) .+ 1
	out = ones(Int, length(ns),2)
	out[1,2] = nd[1]*md
	for i in 2:length(ns)
		out[i,1] = out[i-1,2] + 1
		out[i,2] = out[i,1] + md*nd[i] - 1
	end
	jd = Int((2*S+1)*(2*J+1))
	return ns, nd, out, jd
end
function Hsr(pr,J,S)
	ns, nd, ni, jd = srprep(J,S)
	out = spzeros(Float64,jd,jd)
	out[1:nd[1],1:nd[1]] = Hrot(pr,ns[1]) + Hspi0N(pr,J,S,ns[1])
	for i in 2:length(ns)
		n = ns[i]
		n1part = Hspi1N(pr,J,S,n-1.0)
		out[ni[i-1,1]:ni[i-1,2],    ni[i,1]:ni[i,2]] = n1part
		out[    ni[i,1]:ni[i,2],    ni[i,1]:ni[i,2]] = Hrot(pr,n) + Hspi0N(pr,J,S,n)
		out[    ni[i,1]:ni[i,2],ni[i-1,1]:ni[i-1,2]] = transpose(n1part)
	end
	return out
end
function Htsr(pr,J,S,mcalc,σ)
	md = 2*mcalc + 1
	ns, nd, ni, jd = srprep(J,S,md)
	ni = ni
	ni[1,1] = 1
	out = spzeros(Float64,md*jd,md*jd)
	out[1:ni[1,2],1:ni[1,2]] = kron(eye(md),Hrot(pr,ns[1]) + Hspi0N(pr,J,S,ns[1]))
	out[1:ni[1,2],1:ni[1,2]] += Htor(pr,mcalc,ns[1],σ)
	for i in 2:length(ns)
	n = ns[i]
	n1part = Hspi1N(pr,J,S,n-1.0)
	@inbounds out[ni[i-1,1]:ni[i-1,2],    ni[i,1]:ni[i,2]] = kron(eye(md),n1part)
	@inbounds out[    ni[i,1]:ni[i,2],    ni[i,1]:ni[i,2]] = kron(eye(md),
					Hrot(pr,n) + Hspi0N(pr,J,S,n)) + Htor(pr,mcalc,n,σ)
	@inbounds out[    ni[i,1]:ni[i,2],ni[i-1,1]:ni[i-1,2]] = kron(eye(md),transpose(n1part))
	end
	return out
end

function assign(nmax,j,s,σ,mcalc,mmax,vals,vecs)
	#determine which indices we want to assign
	nlist = Δlist(j,s)
	ndegns = @. 2*nlist + 1
	mcind = mcalc+1
	mcd = 2*mcalc+1
	offset = zeros(Int64,length(nlist))
	offset[2:end] = (mcd) .* ndegns[1:end-1]
	goals = zeros(Int,0,1)
	off = 0
	for n in nlist
        #this will only work for A states
		start = off + Int((mcalc - mmax)*(2*n+1)) + 1
		finsh = off + Int((mcalc + mmax + 1)*(2*n+1))
		part = collect(start:finsh)
		goals = vcat(goals,part)
		off += Int((2*mcalc+1)*(2*n+1))
	end
	assigned = leadfact(vals,copy(vecs))
	vals = vals[assigned]
	vecs = vecs[:,assigned]
	vals = vals[goals]
	vecs = vecs[:,goals]
	#vecs = vecpadder(nlist,ndegns,offset,nmax,mcind,mcd-1,vecs)
	QNs = qngen2(j,s,mmax,σ)
	return QNs, vals, vecs
end

function leadfact(vals,vecs)
	c = zeros(Int,length(vals),2)
	for i in 1:length(vals)
		t = BLAS.iamax(vecs[:,i])
		vecs[t,:] .= 0.0
		c[i,1] = t
		c[i,2] = i
	end
	perm = sortperm(c[:,1])
	return perm
end
function vecpadder(ns,degns,offst,nm,vtmi,vtc,vecs)
	partial = Array{Float64}(undef,0,sum(degns))
	for i in 1:length(ns)
		pad = (nm - ns[i])
		zpad = zeros(Float64, pad, sum(degns))
		for v in 0:(vtc)
			temp = vecs[offst[i]+v*degns[i]+1:offst[i]+(v+1)*degns[i],:]
			temp = vcat(zpad,temp,zpad)
			partial = vcat(partial,temp)
		end
	end
	return partial
end
function jinds(j,s)
	#only seems to work for s=1/2
	snd = convert(Int, s^2 + j*(j+1.0) + (j-s)^2 )
	fnd = convert(Int, s^2 + j*(j+1.0) + (j+s)^2 + 2.0*(j+s) )
	return snd,fnd
end

function tsrdiag(pr,j,s,nmax,mcalc,mmax,σ)
	U = ur(j,s,mcalc)
	if σ==0
		U *= ut(mcalc,j,s)
	end
	H = Matrix(U*Htsr(pr,j,s,mcalc,σ)*U)
	vals, vecs = LAPACK.syev!('V', 'U', H)
	qns, vals, vecs = assign(nmax,j,s,σ,mcalc,mmax,vals,vecs)
	return qns, vals, vecs
end
function tsrcalc(s,nmax,mcalc,mmax,σ)
	prm = paramshift(inprms)
	if isodd(Int64(2*s+1))
		jmin = 0.0
	else
		jmin = 0.5
	end
	jmax = nmax - s
	jarray = collect(Float64,jmin:jmax)
	jd = Int((2*s+1)*sum(2 .* jarray .+ 1))
	outvals = zeros(Float64,Int(jd*(2*mmax+1)))
	outvecs = zeros(Float64,Int((2*s+1)*(2*jmax+2)*(2*mcalc+1)),Int(jd*(2*mmax+1)))
	outqns = zeros(Float64,Int(jd*(2*mmax+1)),6)
	for i in 1:length(jarray)
		j = jarray[i]
		sind, find = jinds(j,s)
		tqns, tvals, tvecs = tsrdiag(prm,j,s,nmax,mcalc,mmax,σ)
		outvals[sind:find] = tvals
		outvecs[1:Int((2*s+1)*(2*j+1)*(2*mcalc+1)),sind:find] = tvecs
		outqns[sind:find,:] = tqns
	end
	return outqns, outvals, outvecs
end

################################################################################
############                        tracalc                         ############
################################################################################

function Tμ(q)
	q = Int(q)
	if q==-1
		return μb/√(2)
	elseif q==0
		return μa
	elseif q==1
		return -μb/√(2)
	else
		return 0.0
	end
end
function intelem(jb,nb,kb,s,j,n,k)
	q = kb-k
	out = √((2*jb+1)*(2*j+1))*WIGXJPF.wig6j(n,j,s,jb,nb,1)*(-1)^(n+s+jb+1)
	out *= √((2*nb+1)*(2*n+1))*WIGXJPF.wig3j(n,1,nb,k,q,-kb)*Tμ(q)*(-1)^(n-kb-1)
	return out
end
function intmat2(jb,jk,s,mcalc,σ,lenb,lenk)
	mat = zeros(Float64,lenk,lenb)
	qnb = qngen2(jb,s,mcalc,σ)
	qnk = qngen2(jk,s,mcalc,σ)
	for y in 1:lenb
	for x in 1:lenk
		@inbounds mat[x,y] = intelem(jb,qnb[y,2],qnb[y,3],s,jk,qnk[x,2],qnb[x,3])
	end
	end
	Uk = ur(jk,s,mcalc)
	Ub = ur(jb,s,mcalc)
	if σ==0
		Uk *= ut(mcalc,jk,s)
		Ub *= ut(mcalc,jb,s)
	end
	mat = Uk*mat*Ub
	return mat
end

function intmat(jb,nb,jk,nk,s,k)
	sqn = length(k)
	mat = zeros(Float64,sqn,sqn)
	for x in 1:sqn
	for y in x:sqn
		@inbounds mat[x,y] = intelem(jb,nb[y],k[y],s,jk,nk[x],k[x])
	end
	end
	return Symmetric(mat)
end
function intbuild(jmax,mcalc,jb,jk,s)
	ns = Int(jmax + s)
	nbs = Δlist(jb,s)
	nks = Δlist(jk,s)
	nbarray = kron(ones((2*mcalc+1)*(2*ns+1)),nbs[1])
	nkarray = kron(ones((2*mcalc+1)*(2*ns+1)),nks[1])
	karray = kron(ones(2*mcalc+1),collect(-ns[1]:ns[1]))
	for i in 2:length(nbs)
		nbarray = vcat(nbarray,kron(ones((2*mcalc+1)*(2*ns+1)),nbs[i]))
		nkarray = vcat(nkarray,kron(ones((2*mcalc+1)*(2*ns+1)),nks[i]))
		karray = vcat(karray,kron(ones(2*mcalc+1),collect(-ns:ns)))
	end
	mat = intmat(jb,nbarray,jk,nkarray,s,karray)
	return mat
end
function intcalc(jmax,mcalc,s)
	sjmd = (2*s+1)*(2*jmax+1)*(2*mcalc+1)
	jind = convert(Int,jmax+s)
	μmat = zeros(Float64,sjmd,sjmd,jind,jind)
	for x in 1:jind
	for y in x:jind
		μmat[:,:,x,y] = intbuild(jmax,y-s,x-s,s)
		μmat[:,:,y,x] = μmat[:,:,x,y]
	end
	end
	return μmat
end

function tracalc(nmax,s,mcalc,σ,qns,vals,vecs)
	jmax = nmax - s
	ns = Δlist(jmax,s)
	trans = zeros(Float64,0,15)
	karray = kron(ones(2*mcalc+1),collect(-ns[1]:ns[1]))
	for i in 2:length(ns)
		karray = vcat(karray,kron(ones(2*mcalc+1),collect(-ns[i]:ns[i])))
	end
	ojb = 1.5
	ojk = 0.5
	lenb = convert(Int,(2*s+1)*(2*ojb+1)*(2*mcalc+1))
	lenk = convert(Int,(2*s+1)*(2*ojk+1)*(2*mcalc+1))
	μmat = intmat2(ojb,ojk,s,mcalc,σ,lenb,lenk)
	for i in 1:length(vals)
	for j in (i+1):length(vals)
		Δj = abs(qns[j,1] - qns[i,1])
		#Δk = abs(qns[j,3] - qns[i,3])
		if Δj ≤ 1
#		Δka = abs(qns[j,3] - qns[i,3])
#		if Δka ≤ 2
			jb = qns[j,1]
			jk = qns[i,1]
			lenb = convert(Int,(2*s+1)*(2*jb+1)*(2*mcalc+1))
			lenk = convert(Int,(2*s+1)*(2*jk+1)*(2*mcalc+1))
			freq = vals[j] - vals[i]
			if (jb!=ojb)||(jk!=ojk)
				μmat = intmat2(jb,jk,s,mcalc,σ,lenb,lenk)
			end
			ojb = jb
			ojk = jk
			if (freq > 0.0)&&(freq < 40.0E+03)#&&(Δk ≤ 3)
				int = (transpose(vecs[1:lenk,j])*μmat*vecs[1:lenb,i])^2# *exp(vals[i]/(TK*2.083661912e+4))
				if int>1.0E-04
					temp = [freq int vals[i]/c transpose(qns[j,:]) transpose(qns[i,:])]
					trans = vcat(trans,temp)
				end
			elseif (freq < 0.0)&&(freq > -40.0E+03)#&&(Δk ≤ 3)
				int = (transpose(vecs[1:lenk,i])*μmat*vecs[1:lenb,j])^2# *exp(vals[j]/(TK*2.083661912e+4))
				if int>1.0E-04
					temp = [-freq int vals[j]/c transpose(qns[i,:]) transpose(qns[j,:])]
					trans = vcat(trans,temp)
				end
			else
			end#freq if
#		else
#		end#Δka if
		else
			break
		end#Δj if
	end#i for
	end#j for
	return trans
end

################################################################################
############                       westerfit                        ############
################################################################################

function qn2ind(n,ka,kc)
	ind = convert(Int,n*(n+1)+ka*(-1)^(n-ka-kc) +1)
end
function qn2ind(j,s,n,ka,kc)
	#only seems to work for s=1/2
	ind = convert(Int,s^2 + j*(j+1.0) + n*(n+1.0) + ka*(-1)^(n-ka-kc))
	return ind
end

function lineprep(lns)
	#converts the input file into a more code friendly format
	#           1  2  3   4   5  6  7  8  9   10 11 12  13   14
	#input  = [ju nu kau kcu mu σu jl nl kal kcl ml σl freq unc]
	#           1  2  3    4  5   6
	#output = [ju σu indu jl σl indl]
	qunus = lns[:,1:12]
	freqs = lns[:,13]
	uncs = lns[:,end]
	inds = zeros(Int64,size(lns)[1],4)
	inds[:,1] = Int64.(2 .* qunus[:,1])
	inds[:,2] = Int64.(2 .* qunus[:,6])
	inds[:,3] = qn2ind.(qunus[:,1],0.5,qunus[:,2],qunus[:,3],qunus[:,4])
	inds[:,4] = Int64.(2 .* qunus[:,7])
	inds[:,5] = Int64.(2 .* qunus[:,12])
	inds[:,6] = qn2ind.(qunus[:,7],0.5,qunus[:,8],qunus[:,9],qunus[:,10])
	#inds = vcat(inds[:,1:2], inds[:,3:4])
	return inds, freqs, uncs
end

function σcount(nfold)
	out = floor(Int,nfold/2)+1
end
function jlister(inds)
	#finds all the unique N values
	jsσs = vcat(inds[:,1],inds[:,4])
	part = vcat(inds[:,2],inds[:,5])
	jsσs = unique(jsσs, dims=1)
	jsσs = jsσs[sortperm(jsσs[:,2]), :]
	return jsσs
end

function limeigcalc(jlist,inds,rp)
	jmax = maximum(jlist[:,1])
	σs = σcount(NFOLD)
	#println("nmax = $nmax")
	evals = zeros(Float64,Int(jd*(2*mmax+1)),σs)
	evecs = zeros(Float64,Int((2*s+1)*(2*jmax+2)*(2*mcalc+1)),Int(jd*(2*mmax+1)),σs)
	eqns = zeros(Float64,Int(jd*(2*mmax+1)),6,σs)
	for i in 1:size(jlist)[1]
		j = jlist[i,1]
		σ = jlist[i,2]
		sind, find = jinds(j,s)
		tqns, tvals, tvecs = tsrdiag(j,s,nmax,mcalc,mmax,σ)
		evals[sind:find,σ+1] = tvals
		evecs[1:Int((2*s+1)*(2*j+1)*(2*mcalc+1)),sind:find,σ+1] = tvecs
		eqns[sind:find,:,σ+1] = tqns
	end
	return evals, evecs, eqns
end
#take difference and determine rms
function rmscalc(vals,inds,ofreqs)
	cfreqs = zeros(size(ofreqs))
	for i in 1:size(cfreqs)[1]
		cfreqs[i] = vals[inds[i,2]] - vals[inds[i,4]]
	end
	omc = ofreqs .- cfreqs
	rms = sqrt(sum(omc .^ 2)/length(omc))
	return rms, omc
end
#construct jacobian
function anaderiv(j,s,σ,mcalc,vec,rp,rpid)
	rp = zeros(Float64, size(rp))
	rp[rpid] = 1.0
	U = ur(j,s,mcalc)
	if σ==0
		U *= ut(mcalc,j,s)
	end
	H = Matrix(U*Htsr(j,s,mcalc,σ)*U)
	vec = vec[1:Int((2*s+1)*(2*j+1)*(2*mcalc+1))]
	out = transpose(vec)*mat*vec
	return out[1]
end
function fderiv()
end
function aderiv()
end
function bderiv()
end
function tθderiv()
end

################################################################################
############                    test parameters                     ############
################################################################################



c = 29979.2458 #MHz to cm-1 converter
F = 5.0*c
ρ = 0.2
V3 = 200.0*c
A = 3000.0
B = 1500.0
C = 1000.0
Dab = 20.0
Fm = 0.0
V6 = 0.0
V3m = 0.0
FN = 0.0
ρm = 0.0
ρ3 = 0.0
FK = 0.0
Fbc = 0.0
Fab = 0.0
V3N = 0.0
V3K = 0.0
V3ab = 0.0
V3bc = 0.0
ρN = 0.0
ρK = 0.0
ρab = 0.0
ρbN = 0.0

ezz = -52.122943
exx = -16.217283
eyy = 0.646016
exz = -24.247198
ezx = -11.934697
ao = -(ezz+eyy+exx)/3.0
a = -(2.0*ezz-eyy-exx)/6.0
d = -(ezx + exz)*0.5
b = (exx - eyy)*0.5

μa = 1.0
μb = 0.5

#=Htot = kron(I(2*mmax+1),Hrot(N)) + Htor(mmax,N)
Hsym = ur(N,mmax)*ut(mmax,N)*Htot*ut(mmax,N)*ur(N,mmax)
λs = eigen(Hsym)
λs.vectors[abs.(λs.vectors) .< 1E-9] .= 0.0=#

#tsrcalc(s,nmax,mcalc,vtmax,σ)
mc = 10
nm = 20
@time mq, ma, mv = tsrcalc(0.5,nm,mc,0,0)
#println(ma ./c)
println(size(ma))

#@time transitions = tracalc(nm,0.5,mc,0,mq,ma,mv)
