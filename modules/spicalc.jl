"""
This constructs the spin-rotation energies and vectors. It is called SpiCalc
Basis set: Ψ=|JNK⟩

Potential improvements:
    Want to remove the f0 variations
"""

module spicalc
using LinearAlgebra
using HalfIntegers
include("./WIGXJPF.jl")
using .WIGXJPF
#Coefficient calculator
#function theta(J::Float64,N::Array{Float64})::Float64
function f(x::Float64,y::Float64)::Float64
    out = sqrt(x*(x+1.0)-y*(y+1.0))
    return out
end
function g(x::Float64,y::Float64)::Float64
    out = sqrt((x-y)*(x-y-1.0))
    return out
end
function theta(j::Float64,n::Float64,s::Float64)::Float64
    if n==0.0
        out = 0.0
    else
        out = n*(n+1.0) + s*(s+1.0) - j*(j+1.0)
        out = out/(2.0*n*(n+1.0))
    end
    return out
end
function phi(j::Float64,n::Float64,s::Float64)::Float64
    if s==0.5
        out = -1.0/n
    else
        out = (n-j+s)*(n+j+s+1)*(s+j-n+1)*(n+j-s)
        out *= 1.0/((2.0*n-1.0)*(2.0*n+1.0))
        out = -sqrt(out)/n
    end
    return out
end
#DelsK = 0.0
#DelsN = 0.0
#DelsNK = 0.0
#DelsKN = 0.0
#delsN = 0.0
#delsK = 0.0

#Matrix Elements
function Hs0K0N(J::Float64,N::Float64,Nt2::Float64,thet::Float64,K::Array{Float64,1},
    srp)::Array{Float64,1}
#spinrotparams = [ao; a; b; 0.0; d; 0.0; DelsN;DelsK; DelsNK; DelsKN; delsN; delsK]
#                  1; 2; 3;   4; 5;   6;    7 ;   8 ;     9 ;    10 ;   11 ;   12 ]
    out  = @. -0.5*srp[1]*(J*(J+1.0)-Nt2-0.75)/thet + (srp[2]*(3.0*K^2 - Nt2) + srp[8]*K^4 + (srp[9] + srp[10])*Nt2*K^2 + srp[7]*Nt2^2)
#    return out
end
function Hs1K0N(J::Float64,N::Float64,Nt2::Float64,K::Array{Float64,1},d::Float64)::Array{Float64,1}
    out = @. d*(K+0.5)*f(N,K)
end
function Hs2K0N(J::Float64,N::Float64,Nt2::Float64,K::Array{Float64,1},
    srp)::Array{Float64,1}
    out = @. (0.5*srp[2]+srp[11]*Nt2+0.5*srp[12]*(K^2+(K+2)^2))*sqrt((Nt2-K*(K-1))*(Nt2-(K-1)*(K-2)))
end
#Matrix Builder
function Hspi0N(J::Float64,N::Float64,S::Float64,srprms)::Array{Float64,2}
    thet = theta(J,N,S)
    Nt2 = N*(N+1.0)
    karray = collect(Float64,-N:N)
    ondiags = Hs0K0N(J,N,Nt2,thet,karray,srprms)
    of1diag = Hs1K0N(J,N,Nt2,karray[2:end],srprms[5])
    of2diag = Hs2K0N(J,N,Nt2,karray[3:end],srprms)
    Spimat = diagm(0=>ondiags,1=>of1diag,2=>of2diag)
    Spimat = Spimat .* fill(thet, size(Spimat))
    Spimat = Hermitian(Spimat)
    return Spimat
end

#List of Matrix Element Generators
#spinrotparams = [ao; a; b; 0.0; d; 0.0; DelsN;DelsK; DelsNK; DelsKN; delsN; delsK]
#                  1; 2; 3;   4; 5;   6;    7 ;   8 ;     9 ;    10 ;   11 ;   12 ]
function Hsm2K1N(N::Float64,K::Array{Float64,1},srp)::Array{Float64,1}
    #out = @. -0.25*(srp[3]+srp[12]*(K*(N-K)+(K-2)*(N-K+2)))*f(N,K-1)*g(N,-K+1.0)
    out = @. -0.25*srp[3]*f(N+1.0,K-2.0)*g(N+1.0,K-2.0)
end
function Hsm1K1N(N::Float64,K::Array{Float64,1},d::Float64)::Array{Float64,1}
    out = @. 0.25*d*(N+2.0*K)*g(N+1.0,K-1.0)
end
function Hs0K1N(N::Float64,K::Array{Float64,1},srp)::Array{Float64,1}
    #out = @. (1.5*srp[2]*K+0.5*K*(srp[8]*K^2+srp[9]*N^2))*sqrt(N^2 - K^2)
    out = @. 1.5*srp[2]*K*sqrt((N+1.0)^2 - K^2)
end
#These off-diag functions need to be accelerated
function Hs1K1N(N::Float64,K::Array{Float64,1},d::Float64)::Array{Float64,1}
    out = @. 0.25*d*(N-2.0*K)*g(N+1.0,-K-1.0)
end
function Hs2K1N(N::Float64,K::Array{Float64,1},srp)::Array{Float64,1}
    #out = @. 0.25*(srp[3]+srp[12]*(K*(N+K)+(K+2)*(N+K+2)))*f(N,K)*g(N,K+1.0)
    out = @. 0.25*srp[3]*f(N+1.0,K+1.0)*g(N+1.0,-K-1.0)
end
#Matrix Builder
function Hspi1N(J::Float64,N::Float64,S::Float64,srprms)::Array{Float64,2}#this one is now structurally corrected
    jphi = -1.0/(J+0.5)
    #N = convert(Float64,N)
    jphi = phi(J,N,S)
    #This function will construct diagonal matrices from the HsnKmN arrays
    #said arrays will be made not square via vcat of zeros
    #these odd matrices will be added together and then outputed
    #remember back when I thought I had fucked these up? I was right!
    karray = collect(Float64,-N:N)
    Spimat = hcat(zeros(Float64,length(karray)-2),diagm(0=>Hs0K1N(N,karray[2:end-1],srprms)),zeros(Float64,length(karray)-2))
    #The next two lines are very expensive compuationally
    Spimat += hcat(diagm(0=>Hs1K1N(N,karray[1:end-2],srprms[5]),-1=>Hsm2K1N(N,karray[4:end],srprms)),zeros(Float64,length(karray)-2,2))
    Spimat += hcat(zeros(Float64,length(karray)-2,2),diagm(0=>Hsm1K1N(N,karray[3:end],srprms[5]),1=>Hs2K1N(N,karray[1:end-3],srprms)))
    Spimat = Spimat .* fill(jphi, size(Spimat))
    return Spimat
end

#this idea is super gross. Sorry
#function welcomesr(parameters)
#    eaa = parameters[1, 1]
#    ebb = parameters[2, 2]
#    ecc = parameters[3, 3]
#    ezx = parameters[1, 2]
#    exz = parameters[2, 1]
#    ao = -(eaa+ecc+ebb)/3.0
#    a = -(2.0*eaa-ecc-ebb)/6.0
#    d = -(ezx + exz)*0.5
#    b = -(ebb-ecc)*0.5
#    global a
#    global ao
#    global b
#    global d
#    cow = 0.0
#end

ϵzz = 0.0
ϵxx = 0.0
ϵyy = 0.0
ϵxz = 0.0
ϵzx = 0.0

#Tensor version!
ϵtns = zeros(Float64,3,5)
ϵtns[1,3] = -(ϵxx+ϵyy+ϵzz)/sqrt(3.0)
ϵtns[2,2] =  0.5*(ϵxz-ϵzx)#T^1_-1
ϵtns[2,3] =  0.0#T^1_0
ϵtns[2,4] =  0.5*(ϵxz-ϵzx)#T^1_1
ϵtns[3,1] =  0.5*(ϵxx-ϵyy)#T^2_-2
ϵtns[3,2] =  0.5*(ϵxz+ϵzx)#T^2_-1
ϵtns[3,3] =  (2.0*ϵzz-ϵxx-ϵyy)/sqrt(6.0)#T^2_0
ϵtns[3,4] =  0.5*(ϵxz+ϵzx)#T^2_1
ϵtns[3,5] =  0.5*(ϵxx-ϵyy) #T^2_2

function srelems(ϵtns,J,S,N,Ki,q,Nb)
    out = 0.0
    kf = collect(Float64,0:2)
    for k in 0:2
        @. p0  = sqrt(N*(N+1.0)*(2.0*N+1.0))*wig6j(1,1,k,Nb,N,N)
        @. p0 += sqrt(Nb*(Nb+1.0)*(2.0*Nb+1.0))*wig6j(1,1,k,N,Nb,Nb)
        @. p0 *= (2.0*kf[k+1]+1.0)*(-1.0)^(Nb-Ki-q)*wig3j(Nb,ki,N,-Ki-q,q,Ki)
        @. p0 *= ϵtns[k+1,q+3]
        out += p0
    end
    #This line should be moved to the matrix function
    out *= 0.5*sqrt(S*(S+1.0)*(2.0*S+1.0)*(2.0*N+1.0)*(2.0*Nb+1.0))
    out *= wig6j(N,S,J,S,Nb,1)*(-1.0)^(J+S+Nb)
end

function srmat(ϵtns,J,S,RotMats)
    nmin = J - S
    nmax = J + S
    narray = vcat(0,collect(Int64, nmin:nmax))
    nfarray = vcat(0.0,collect(Float64, nmin:nmax))
    mat = zeros(Float64, sum(2 .* narray .+ 1),sum(2 .* narray .+ 1))
    for n in 2:(nmax-min+2)
        N = narray[n]
        N = nfarray[n]
        karray = collect(Int64,-n:n)
        for nb in N:(nmax-min+1)
            Nb = narray[nb]
            Nb = nfarray[nb]
            if n==nb
                sind = sum(2 .* narray[1:n-1] .+ 1)
                eind = sum(2 .* narray[1:n] .+ 1)
                mat[sind:eind,sind:eind] = rotcall(n,RotMats)
                for q in -2:2
                    mat[sind:eind,sind:sind] .+= diagm(q=>srelems(ϵtns,J,S,N,karray[abs(q)+1:end],q,Nb,J,S,N,Nb))
                end
            else
                for q in -2:2
                    mat[sind:eind,(sind+2*nb-1):(end+2*nb-1)] .+= diagm(q=>srelems(ϵtns,J,S,N,karray[abs(q)+1:end],q,Nb,J,S,N,Nb))
                end
            end
        end
    end
    return Symmetric(mat)
end

function buildall(ϵtns,Jmax,S,RotMats)
    if typeof(S)==Int64
        Jcap = Jmax
        Jmin = 0
    else
        Jcap = Int64(Jmax+S)
        Jmin = 1//2
    end
    mat = zeros(Float64,2*Jmax+1,2*Jmax+1,Jcap)
    jfs = collect(Float64,1:Jcap)
    S = Float64(S)
    for j in Jmin:Jcap
        mat[:,:,j] = srmat(ϵtns,j,S,jfs[j],S,RotMats)
    end
    return mat
end

end
