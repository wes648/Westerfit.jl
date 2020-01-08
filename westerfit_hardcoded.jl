"""
This will serve as a test form for TorSpiSymCalc with a hardcoded Hamiltonian
A generalized from of SymCalc currently exists. This should make the creation
of TorSpiSymFit easier. Or at least I hope it will

Some efficiency things to work on: It currently takes abour 4 second to
calculate a Jmax = 4.5 which is far too slow.
-Variable Types
    removing ambiguities should speed up the code
-@. syntax should provide some improvements
    need to read much more on this
-condensable for loops
-@fastmath should improve things as we have no complex values

Also should consider making matrix type specific functions
Hrot(N)
Hspi(J,N)
Hspion1(J,N)
Htor(J,N,vt,sigma)
Htorov1(J,N,vt,sigma)
Htorov2(J,N,vt,sigma)
"""

println("")
using LinearAlgebra
using LinearAlgebra.BLAS
using Profile
using DelimitedFiles
#using Traceur

#Etors = zeros(5,3)
useSpi=true
useTor=false
useNatUnits=false
A = 0.0
B = 0.0
C = 0.0
Dab = 0.0
DelNK = 0.0
DelK = 0.0
DelN = 0.0
delN = 0.0
delK = 0.0
F = 0.0
V3 = 0.0
V6 = 0.1
rho = 0.0
fv3 = 0.0
fv6 = 0.0
ezz = 0.0
exx = 0.0
eyy = 0.0
exz = 0.0
ezx = 0.0
iota = 1.0
μa = 1.0
μb = 0.8
vtmax = 0
vtcalc = 4
Nmax = 0
vthardcap = 10

#You may have noticed the absurd number of νt limits. Here is what they mean in the order they are used in the code:
#νthardcap is the value used in the 1st diagonalization stage. BELGI uses 10. It has a silly name because it's probably not worth futzing with unless you are calculating very high torsional states
#vtcalc is the highest vt in the 2nd diagonalization stage. The truncation is the improve calculation speeds. Should probably be at least 3 higher than what you want in the output
#vtmax is the highest vt that is output by the program.

print("Enter molecule name: ")
#molnam = readline()
molnam = "test"
include("$molnam.inp")
ρ=rho
useSpi=true
#const Jmax = Jamx
#const vtmax = vtmax
#const vtcalc = vtcalc
if useNatUnits==true
    F = 29979.2458*F
    V3 = 29979.2458*V3
else
    nothing
end

if useSpi==true
    ao = -(ezz+eyy+exx)/3.0
    a = -(2.0*ezz-eyy-exx)/6.0
    d = -(ezx + exz)*0.5
    b = (exx-eyy)*0.5
else
    ao = 0.0
    a = 0.0
    d = 0.0
    b = 0.0
end
if useTor==false
    vtmax = 0
    vtcalc = 0
    ρ = 0.0
    F = 0.0
    V3 = 0.0
else
    thetaRAM = 0.5*atan(2*Dab/(A-B))
    μa = μa*cos(thetaRAM) + μb*sin(thetaRAM)
    μb = -μa*sin(thetaRAM) + μb*cos(thetaRAM)
end

## Below are some general purpose functions ##
function f(x,y)::Float64
    out = sqrt(x*(x+1.0)-y*(y+1.0))
    return out
end
function g(x,y)::Float64
    out = sqrt((x-y)*(x-y-1))
    return out
end
function theta(J::Float64,Nf)::Float64
    out = J*(J+1.0) - Nf*(Nf+1.0) - 0.75
    out = -out/(2.0*Nf*(Nf+1.0))
    return out
end
function vecvalpair(m)
#This function provides the definitional pairing of the torsional eigenvalues to the eigenvector elements
	out = @. 11 - div((m-(m%2))*((-1)^(m%2)),2)
#    out = @. 11 + out
    return out
end


########################################################
### Here is the purely torsional protion of H_{tot} ###
#######################################################

function Htorm0(K,m,sigma)
    out = @. F*(3.0*m+sigma-rho*K)^2 + V3*0.5 + V6*0.5 + 0.0*fv3*(3.0*m+sigma)^2 + 0.0*fv6*(3.0*m+sigma)^2
end
function Htorm1(K,m,sigma)
    out = @. -V3*0.25 + 0.0*fv3*(3.0*m+sigma)^2
end
function Htorm2(K,m,sigma)
    out = @. -V6*0.25 + 0.0*fv6*(3.0*m+sigma)^2
end
function Htors(K,sigma)
    marray = collect(-vthardcap:vthardcap)
    ondiags = Htorm0(K,marray,sigma)
    of1diag = Htorm1(K,marray[1:end-1],sigma)
    of2diag = Htorm2(K,marray[1:end-2],sigma)
    Tormat = diagm(0=>ondiags, 1=>of1diag, 2=>of2diag)
    Tormat = Symmetric(Tormat)
    Tormat = eigen!(Tormat)
#   #println(Tormat.values)
#   return Tormat.values, Tormat.vectors
    return Tormat
end
function tormatch(eigs)#vals,vecs)#,p)
#    println(p)
#    matched = @. findmax(eigs.vectors[:,p])[2]
#    println(matched)
#   #matched = zeros(Int64,2*vthardcap+1)
#   #second = zeros(Int64,2*vthardcap+1)
    assigned = zeros(Float64,2*vthardcap+1)
    reordered = zeros(Float64,(2*vthardcap+1,2*vthardcap+1))
#   #parray = collect(Int64,1:2*vthardcap)
    test = collect(Int8,1:2*vthardcap+1)
    newindices = vecvalpair(test)
    eigs.values[newindices] = eigs.values[test]
    eigs.vectors[:,newindices] = eigs.vectors[:,test]
    return eigs#assigned,reordered
end
function TorCalc(Nm,sigma)
    #Okay so I need to reorder the eigenvectors as well.
    #I will out put a 3D array of eigenvectors and a 2D of values.
    #This will not be truncated. Truncation will occur in what is called in the next stage
    #marray = collect(1:21)
    toreigs = Htors(-Nm,sigma)
    toreigs = tormatch(toreigs)#tores,torvs)#,marray)
    torvals = toreigs.values
    torvecs = toreigs.vectors
    for k in (-Nm+1):Nm
        neweigs = Htors(k,sigma)
        neweigs = tormatch(toreigs)#tores,torvs)#,marray)
        torvals = hcat(torvals,neweigs.values)
        torvecs = cat(torvecs,neweigs.vectors,dims=3)
    end
    return torvals,torvecs
end


#################################################
### Here is the rotational protion of H_{tot} ###
#################################################
function Hr0K(N,Nt2,K)
    out = @. (0.5*(B+C)+DelNK*K^2+DelN*Nt2)*Nt2+(A-0.5*(B+C)+DelK*K^2)*K^2
end
function Hr1K(N,Nt2,K)
    out = @. Dab*sqrt(Nt2-K*(K+1))*(K+0.5)
end
function Hr2K(N,Nt2,K)
    out = @. sqrt((Nt2-K*(K+1))*(Nt2-(K+1)*(K+2)))*(0.25*(B-C) + delN*Nt2 + 0.5*delK*(K^2 + (K+2)^2))
end
function Hrot(N)
    Nt2 = Float64(N*(N+1.0))
    karray = collect(-N:N)
    ondiags = Hr0K(N,Nt2,karray)
    of1diag = Hr1K(N,Nt2,karray[1:end-1])
    of2diag = Hr2K(N,Nt2,karray[1:end-2])
    Rotmat = diagm(0=>ondiags,1=>of1diag,2=>of2diag)
    Rotmat = Symmetric(Rotmat)
    return Rotmat
end

#####################################################
### Here is the spin-rotation protion of H_{tot} ###
####################################################

function Hs0K0N(J,N,Nt2,thet,K)
    out  = @. -0.5*ao*(J*(J+1.0)-Nt2-0.75) + a*(3.0*K^2 - Nt2)*thet
#    return out
end
function Hs1K0N(J,N,Nt2,thet,K)
    out = @. d*(K+0.5)*f(N,K)*thet
end
function Hs2K0N(J,N,Nt2,thet,K)
    out = @. 0.5*b*f(N,K)*f(N,K+1)*thet
end
function Hspi0N(J,N)
    thet = theta(J,N)
    Nt2 = Float64(N*(N+1.0))
    karray = collect(-N:N)
    ondiags = Hs0K0N(J,N,Nt2,thet,karray)
    of1diag = Hs1K0N(J,N,Nt2,thet,karray[1:end-1])
    of2diag = Hs2K0N(J,N,Nt2,thet,karray[1:end-2])
    Spimat = diagm(0=>ondiags,1=>of1diag,2=>of2diag)
    Spimat = Symmetric(Spimat)
    return Spimat
end
function Hspi0Nf0(J,N)
    Nt2 = Float64(N*(N+1.0))
    karray = collect(-N:N)
    ondiags = Hs0K0N(J,N,Nt2,0.5,karray)
    of1diag = Hs1K0N(J,N,Nt2,0.5,karray[1:end-1])
    of2diag = Hs2K0N(J,N,Nt2,0.5,karray[1:end-2])
    Spimat = diagm(0=>ondiags,1=>of1diag,2=>of2diag)
    Spimat = Symmetric(Spimat)
    return Spimat
end
function Hs0K1N(jphi,N,K)
    out = @. 1.5*a*K*sqrt(N^2 - K^2)*jphi
end
function Hs1K1N(jphi,N,K)
    out = @. 0.25*d*(N+2.0*K+1.0)*g(N,K)*jphi
end
function Hsm1K1N(jphi,N,K)
    out = @. 0.25*d*(N-2.0*K+1.0)*g(N,-K)*jphi
end
function Hs2K1N(jphi,N,K)
    out = @. 0.25*b*f(N,K)*g(N,K+1.0)*jphi
end
function Hsm2K1N(jphi,N,K)
    out = @. -0.25*b*f(N,K-1)*g(N,-K+1.0)*jphi
end
function Hspi1N(J,N)
    jphi = -sqrt((2*J+2)/(2*J+1))/(J+0.5)
    Nf = Float64(N)
    #This function will construct diagonal matrices from the HsnKmN arrays
    #said arrays will be made not square via vcat of zeros
    #these odd matrices will be added together and then outputed
    karray = collect(-(N):(N))
    Spimat = vcat(zeros(1,length(karray)-2),diagm(0=>Hs0K1N(jphi,Nf,karray[2:end-1])),zeros(1,length(karray)-2))
    Spimat += vcat(diagm(0=>Hs1K1N(jphi,Nf,karray[1:end-2])),zeros(2,length(karray)-2))
    Spimat += vcat(zeros(2,length(karray)-2),diagm(0=>Hsm1K1N(jphi,Nf,karray[3:end])))
    #println(size(Spimat))
    p2Ks = diagm(1=>Hs2K1N(jphi,Nf,karray[1:end-3]))
    #println(size(p2Ks)[1])
    #println(length(karray))
    Spimat += vcat(p2Ks,zeros(2,length(karray)-2))
    m2Ks = diagm(-1=>Hsm2K1N(jphi,Nf,karray[4:end]))
    #println(size(m2Ks)[1])
    Spimat += vcat(zeros(2,length(karray)-2),m2Ks)
    return Spimat
end
function Hspi1Nf0(J,N)
    jphi = -sqrt((2*J+2)/(2*J+1))/(J+0.5)
    Nf = N
    karray = collect(-(N):(N))
    Spimat = vcat(zeros(1,length(karray)-2),diagm(0=>Hs0K1N(jphi,Nf,karray[2:end-1])),zeros(1,length(karray)-2))
    Spimat += vcat(diagm(0=>Hs1K1N(jphi,Nf,karray[1:end-2])),zeros(2,length(karray)-2))
    Spimat += vcat(zeros(2,length(karray)-2),diagm(0=>Hsm1K1N(jphi,Nf,karray[3:end])))
    return Spimat
end

#function HspiRotN(J,N,Nlmat=0.0)
#    Numat = Hrot(N)+Hspi0N(J,N)
#    URmat = Hspi1N(J,N)
#    llmat = zeros(2*N-1,2*N+1)
#    Hsrmat = [Numat URmat; llmat Nlmat]
#    Hsrmat = Symmetric(Hsrmat)
#    return Hsrmat, Numat
#end


########################################################
### Here is the coupled-torsional protion of H_{tot} ###
########################################################

function Ht0K0m0N(tvals,tvecs,thet,J,N,K,vt,sigma)
    out = zeros(Float64,2*vthardcap+1)
    for m in -vthardcap:vthardcap
        mp = m+vthardcap+1
        out[mp] = tvals[mp,K+Nmax+1]
        out[mp] += iota*(3.0*m+sigma)*K*thet
        out[mp] *= tvecs[mp,vt+1,K+Nmax+1]^2
    end
    out = sum(out)
    return out
end
function Ht0K0m1N(jphi,J,N,K,vt,sigma,tvecs)
    out = zeros(Float64,2*vthardcap+1)
    for m in -vthardcap:vthardcap
        mp = m+vthardcap+1
        out[mp] = iota*(3.0*m+sigma)
        out[mp] *= sqrt(N^2 - K^2)*jphi
        out[mp] *= tvecs[mp,vt+1,K+Nmax+1]*tvecs[mp,vt+2,K+Nmax+1]
    end
    out = sum(out)
    return out
end

#The below functions are only for higher order terms that I can turn on later
#function Ht1K0m0N(J,N,vt,sigma)
#function Ht2K0m0N(J,N,vt,sigma)
#function Ht0K1m0N(J,N,vt,sigma)
#function Ht1K1m0N(J,N,vt,sigma)
#function Ht2K1m0N(J,N,vt,sigma)
#function Ht0K2m0N(J,N,vt,sigma)
#function Ht1K2m0N(J,N,vt,sigma)
#function Ht2K2m0N(J,N,vt,sigma)
#function Ht1K0m1N(J,N,vt,sigma)
#function Ht2K0m1N(J,N,vt,sigma)
#function Ht0K1m1N(J,N,vt,sigma)
#function Ht1K1m1N(J,N,vt,sigma)
#function Ht2K1m1N(J,N,vt,sigma)
#function Ht0K2m1N(J,N,vt,sigma)
#function Ht1K2m1N(J,N,vt,sigma)
#function Ht2K2m1N(J,N,vt,sigma)

function Hspitor0N(J,N,v,sigma,tvals,tvecs)
    thet = theta(J,N)
    #karray = collect(Int64,-N:N)
    elems = zeros(Float64,2*N+1)
    for k in -N:N
        elems[k+N+1] = Ht0K0m0N(tvals,tvecs,thet,J,N,k,v,sigma)
    end
    out = diagm(0=>elems)
    return out
end
function Hspitor1N(J,N,m,sigma,tvecs)
    jphi = -sqrt((2*J+2)/(2*J+1))/(J+0.5)
    Nf = Float64(N)
    karray = collect(-(N):(N))
    elems = zeros(Float64,2*N-1)
    for k in (-N+1):(N-1)
        elems[k+N] = Ht0K0m1N(jphi,J,Nf,k,m,sigma,tvecs)
    end
    tsmat = vcat(zeros(1,length(karray)-2),diagm(0=>elems),zeros(1,length(karray)-2))
    return tsmat
end
function Htot(J,N,sigma,tvals,tvecs)
    thet = theta(J,N)
    marray = collect(-vthardcap:vthardcap)
    Hsrumat = Hrot(N) + Hspi0N(J,N)
    Numat = Hsrumat + Hspitor0N(J,N,0,sigma,tvals,tvecs)
    Hsrlmat = Hrot(N-1) + Hspi0N(J,N-1)
    Nlmat = Hsrlmat + Hspitor0N(J,N-1,0,sigma,tvals,tvecs)
    URsrmat = Hspi1N(J,N)
    URmat = URsrmat + Hspitor1N(J,N,0,sigma,tvecs)
    for vt in 1:vtcalc
        #I should figure out how to vectorize this but that's a problem for future Wes
        vtumat = Hsrumat + Hspitor0N(J,N,vt,sigma,tvals,tvecs)
        Numat = cat(Numat, vtumat,dims=(1,2))
        vtlmat = Hsrlmat + Hspitor0N(J,N-1,vt,sigma,tvals,tvecs)
        Nlmat = cat(Nlmat,vtlmat;dims=(1,2))
        vturmat = URsrmat + Hspitor1N(J,N,vt,sigma,tvecs)
        URmat = cat(URmat,vturmat;dims=(1,2))
    end
    LLmat = transpose(zeros(size(URmat)))
    Htotmat = [Numat URmat; LLmat Nlmat]
    Htotmat = Symmetric(Htotmat)
    return Htotmat
end
#fuck I need to fix this
function Htotf0(J,sigma,tvals,tvecs)
    thet = theta(J,1)
    Hsrumat = Hrot(1) + Hspi0Nf0(J,1)
    Numat = Hsrumat + Hspitor0N(J,1,0,sigma,tvals,tvecs)
    #Hsrlmat = 0.0 #Hrot(N-1) + Hspi0N(J,N-1)
    Nlmat = 0.0 #Hsrlmat + Hspitor0N(J,N-1,0,sigma)
    URsrmat = Hspi1Nf0(J,1)
    URmat = URsrmat + Hspitor1N(J,1,0,sigma,tvecs)
    for vt in 1:vtcalc
        #I should figure out how to vectorize this but that's a problem for future Wes
        vtumat = Hsrumat + Hspitor0N(J,1,vt,sigma,tvals,tvecs)
        Numat = cat(Numat, vtumat,dims=(1,2))
        vtlmat = 0.0#Hsrlmat + Hspitor0N(J,N-1,vt,sigma)
        Nlmat = cat(Nlmat,vtlmat;dims=(1,2))
        vturmat = URsrmat + Hspitor1N(J,1,vt,sigma,tvecs)
        URmat = cat(URmat,vturmat;dims=(1,2))
    end
    LLmat = transpose(zeros(size(URmat)))
    Htotmat = [Numat URmat; LLmat Nlmat]
    Htotmat = Symmetric(Htotmat)
    return Htotmat
end

#now i need a function to diagonalize and assign quantum numbers
#find index of highest value in eigenvector, assign QNs as function of index
#I'm going to just use bins and then definitionally assign Ks
#the bins will be temorary 2D arrays, one for each N
#then each column will share an vt value
#finally each column will be sorted to assign the K
function qnassign(J,sigma,eigs)
    Nu = Int(J+0.5)
    Nl = Int(J-0.5)
    vtcalcd = vtcalc+1
    Nueigs = zeros(2*Nu+1,vtcalcd)
    Nleigs = zeros(2*Nl+1,vtcalcd)
    for i in 1:length(eigs.values)
        energy = eigs.values[i]
        #println(energy)
        #we are using the highest value in the eigenvector to understand where in
        #the original matrix each eigenvalue came from to figure out what it's QNs must be
        index = iamax(eigs.vectors[:,i])
        #this if statement determines if N = J±1/2
        if index<=((vtcalcd)*(2*J+2))
            N = Nu
            #Then we use it to determine vt
            vtp = convert(Int,floor(index/(2*N+1.01)))+1
            Nueigs[findfirst(isequal(0.0),Nueigs[:,vtp]),vtp] = energy
        else
            N = Nl
            index = index-(vtcalcd)*(2*J+2)
            vtp = convert(Int,floor(index/(2*N+1.01)))+1
            Nleigs[findfirst(isequal(0.0),Nleigs[:,vtp]),vtp] = energy
        end
    end
    #sort by energy to "assign" K
    vmaxind = vtmax+1
#    vmaxind = 2*(vtcalc-vtmax)
#    vminind = 2*(vtcalc-2*vtmax)
#    println("$vmaxind, $vminind")
#    println(size(Nueigs))
#    println(Nueigs)
    Nueigs = sort(Nueigs[:,1:vmaxind],dims=1)
    Nleigs = sort(Nleigs[:,1:vmaxind],dims=1)
    Jeigs = vcat(Nleigs[:,1:vmaxind],Nueigs[:,1:vmaxind])
    return Jeigs
end

function TSRDiag(N,sigma,tvals,tvecs)
    J = N-0.5
    Hmat = Htot(J,N,sigma,tvals,tvecs)
    Hmat = eigen!(Hmat)
    output = qnassign(J, sigma, Hmat)
    return output
end

function TSRDiagf0(sigma,tvals,tvecs)
    J = 0.5
    Hmat = Htotf0(J,sigma,tvals,tvecs)
    Hmat = eigen!(Hmat)
    output = qnassign(J, sigma, Hmat)
    return output
end

function RotSpiTorCalc(Nmax,sigma,tvals,tvecs)
    marray = [-vthardcap:vthardcap]
    energies = TSRDiagf0(sigma,tvals,tvecs)
    for n = 2:Nmax
        #println("n=$n")
        news = TSRDiag(n,sigma,tvals,tvecs)
        energies = vcat(energies,news)
    end
    return energies
end

#Ators = TorCalc(Nmax+1,0)
#E1tors = TorCalc(Nmax+1,1)
#E2tors = TorCalc(Nmax+1,-1)
#TorEs = cat(E2tors,Ators,E1tors,dims=3)
#println(size(TorEs))
#U = UniformScaling(2)
#Etors = U \ (E1tors+E2tors)
#@time a=RotSpiTorCalc(Nmax,0)
#println(size(a))
#@time b=RotSpiTorCalc(Nmax,-1)
#println(size(b))
#@time c=RotSpiTorCalc(Nmax,1)
#println(size(c))

function fullandfinal()
    σ0val,σ0vec = TorCalc(Nmax+1,0)
    σ1val,σ1vec = TorCalc(Nmax+1,1)
    #E2tors = TorCalc(Nmax+1,-1)
    #alltors = cat(E2tors,Ators,E1tors,dims=3)
    sig0ful=RotSpiTorCalc(Nmax,0,σ0val,σ0vec)
#    println(size(sig0ful))
    sig1ful=RotSpiTorCalc(Nmax,1,σ1val,σ1vec)
#    println(size(sig1ful))
    #sigm1ful=RotSpiTorCalc(Nmax,-1)
    #println(size(sigm1ful))
    finaleigs = cat(sig0ful,sig1ful,dims=3)
end
#There is still plenty of optimizations to do here
#also need to figure out what's up with σ=±1
"""
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



function TraCalc(Nmax,energies,sigma,vt)
    #much of the following array creation should be moved into a separate function as TraCalc
    #This is because it will be the same for each vt and σ pair. Doing this first could open up parallel options
    arr = collect(1:Nmax-1)
    arr = collect(Iterators.flatten(zip(arr,arr)))
    #this excluded list is for the levels that won't have a-type transitions
    exclude = zeros(Int64,length(arr))
    exclude[1] = 2
    for i=1:(length(arr)-1)
        exclude[i+1] = exclude[i]+2*arr[i]+1
    end
#    println(exclude)
    #this list has the J value for the level of identical index
    posjs = collect(0.5:Nmax-0.5)
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
#    println(ninds)
#    println(size(jinds))
#    println(size(ninds))
#   freqs is actually frequencies[1] and intensities[2] and type[3]
    freqs = Array{Float64}(undef,0,3)
#   states contains the indices of the states for a given frequency
    states = Array{Int64}(undef,0,2)
    leng = size(energies)[1]
#    println(length(posjs))
#    println(length(ninds))
#    println(length(kinds))
#    println(length(jinds))
#    println(size(energies))
#    println(leng)
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
            intent = (N+1+K)*(N+1-K)/(N+1)
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
            intent = (2*N+1)*K^2/(N^2+N)
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
            intent = (N+1+K)*(N+1-K)/(N+1)
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
        if i in exclude
            nothing
        else
        #   calcas
            #ΔJ= 0,ΔN=+1,ΔK=-1  a
            indexu = i+2*J
            indexu = convert(Int64,indexu)
            if (indexu≤leng)&&(J-N>0)
                freq = energies[indexu,vtind,sigind]-energies[i,vtind,sigind]
                intent = (N+2-K)*(N+1-K)/(4*N+4)
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
                intent = (N+1-K)*(N+K)*(2*N+1)/(4*N*(N+1))
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
                intent = (N+2-K)*(N+1-K)/(4*N+4)
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
    end
    #At this point there are two vital arrays: freqs and states
    #println(sigma==0)
    #println(size(states))
    #println(size(freqs))
    #println(states)
    #println(ninds)
    #println(ninds[states[:,1]])
    if sigma==0
        ju = jinds[states[:,1]]
        jl = jinds[states[:,2]]
        nu = ninds[states[:,1]]
        nl = ninds[states[:,2]]
        ku = kinds[states[:,1]]
        kl = kinds[states[:,2]]
        sigs = fill(sigma, length(ju))
        open("./AStates.txt", "w") do io
            writedlm(io, [ju nu ku sigs jl nl kl sigs freqs])
        end
    else
        ju = jinds[states[:,1]]
        jl = jinds[states[:,2]]
        nu = ninds[states[:,1]]
        nl = ninds[states[:,2]]
        ku = kinds[states[:,1]]
        kl = kinds[states[:,2]]
        sigs = fill(sigma, length(ju))
        open("EStates.txt", "w") do io
            writedlm(io, [ju nu ku sigs jl nl kl sigs freqs])
        end
    end
    return freqs, states
end


#Nmax = 4
#println(Nmax)

@time finalenergies=fullandfinal()

@time Afrequencies, Aqnindices = TraCalc(Nmax,finalenergies,0,0)
@time Efrequencies, Eqnindices = TraCalc(Nmax,finalenergies,1,0)


println()
println("A mircle has come to pass. This code seems to have properly run")
println()
