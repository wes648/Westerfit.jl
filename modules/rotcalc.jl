"""
This constructs the rotational energy matrices. It is called RotCalc
Basis set: Ψ=|NK⟩

Potential improvements:
    complete dictionary conversion of parameters
    not sure of specifics but I'm sure of plenty
"""

module rotcalc
using LinearAlgebra
include("./WIGXJPF.jl")
using .WIGXJPF
#once again questioning the signs in my N_± elements
#function Hr0K(N::Float64,Nt2::Float64,K::Array{Float64,1},params)::Array{Float64,1}
#rotparams = [A; B; C; Dab; DelN; DelNK; DelK; delN; delK]
#             1; 2; 3;  4;   5  ;   6  ;   7;    8 ;  9
function Hr0K(N::Float64,Nt2::Float64,K::Array{Float64,1},rprms)::Array{Float64,1}
    out = @. (0.5*(rprms[2] + rprms[3]) + rprms[6]*K^2 + rprms[5]*Nt2)*Nt2
    out += @. (rprms[1]-0.5*(rprms[2] + rprms[3]) + rprms[7]*K^2)*K^2
end
#function Hr1K(N::Float64,Nt2::Float64,K::Array{Float64,1},params)::Array{Float64,1}
function Hr1K(N::Float64,Nt2::Float64,K::Array{Float64,1},Dab::Float64)::Array{Float64,1}
    out = @. Dab*sqrt(Nt2-K*(K+1))*(K+0.5)
end
#function Hr2K(N::Float64,Nt2::Float64,K::Array{Float64,1},params)::Array{Float64,1}
function Hr2K(N::Float64,Nt2::Float64,K::Array{Float64,1},rprms)::Array{Float64,1}
    out = @. 0.25*(rprms[2] - rprms[3]) + rprms[8]*Nt2 + 0.5*rprms[9]*(K^2 + (K-2)^2)
    out = @. out*sqrt((Nt2 - K*(K - 1))*(Nt2 - (K - 1)*(K - 2)))
end

function Hrot(Nf::Float64,rotparams)
    Nt2 = Nf*(Nf+1.0)
#    N = convert(Float64,N)
    karray = collect(Float64,-Nf:Nf)
#    ondiags = Hr0K(Nf,Nt2,karray,parameters)
#    of1diag = Hr1K(Nf,Nt2,karray[2:end],parameters)
#    of2diag = Hr2K(Nf,Nt2,karray[3:end],parameters)
    ondiags = Hr0K(Nf,Nt2,karray,rotparams)
    of1diag = Hr1K(Nf,Nt2,karray[2:end],rotparams[4])
    of2diag = Hr2K(Nf,Nt2,karray[3:end],rotparams)
    Rotmat = diagm(0=>ondiags,1=>of1diag,2=>of2diag)
    Rotmat = Hermitian(Rotmat)
    return Rotmat
end

function RotBuild(parameters, Nmax)
    #global A = parameters[1, 1]
    #global B = parameters[1, 2]
    #global C = parameters[1, 3]
    #global DelN = parameters[2, 1]
    #global DelK = parameters[2, 2]
    #global DelNK = parameters[2, 3]
    #global delN = parameters[3, 1]
    #global delK = parameters[3, 2]
    #global Dab = parameters[3, 3]
    bigmat = zeros(Float64, 3, 3)
    counter = 3
    for n in 1:Nmax
        rotmat = Hrot(parameters,n)
        bigmat = cat(bigmat, rotmat, dims=3)
        bigmat = hcat(bigmat, zeros(Float64, counter, 2, n+1))
        counter += 2
        bigmat = vcat(bigmat, zeros(Float64, 2, counter, n+1))
    end
end

#=
#Okay here are some tensor versions
Btns = zeros(Float64, 3,5) #B^k_q+3
Btns[1,3] = -(A+B+C)/sqrt(3.0)#B^0_0
Btns[3,5] = 0.5*(B-C)#B^2_2
Btns[3,4] = -0.5*Dab#B^2_1
Btns[3,3] = (2*A-B-C)/sqrt(6.0)#B^2_0
Btns[3,2] = -0.5*Dab#B^2_-2
Btns[3,1] = 0.5*(B-C)#B^2_-2
=#

function rote(Btns,N,K,q)
    #K is the ket value
    #q is the difference in the ket values
    #-Kb + K + q = 0; Kb = K + q
    #k=0
    out = @. Btns[1,3+q]*wig3j(N,0,N,-K-q,q,K)*wig6j(N,N,1,0,1,N)
    #k=2
    out += @. Btns[3,3+q]*wig3j(N,2,N,-K-q,q,K)*sqrt(5.0)*wig6j(N,N,1,2,1,N)
    out = @. out .* (-1.0)^(N-K-q) #Reworking this line could have minor speed boosts
    return out
end

function rotmat(Btns,Ni,Nf)
    karray=collect(Int64,-Ni:Ni)
    mat = zeros(Float64, 2*Ni+1,2*Ni+1)
    for q in -2:2
        mat += diagm(q=>rote(Btns,Ni,karray[abs(q)+1:end],q))
    end
    mat = mat .* Nf*(Nf+1.0)*(2.0*Nf+1.0)
    return mat
end

function rotcall(Ni::Int64,bigmat::Array{Float64,3})::Array{Float64,2}
    if Ni==0
        return zeros(Float64, 1,1)
    else
        return bigmat[1:(2*Ni),1:(2*Ni),Ni]
    end
end

function buildall(Btns,Nmax)
    mat = zeros(Float64,2*Nmax+1,2*Nmax+1,Nmax)
    nfs = collect(Float64,1:Nmax)
    for n in 1:Nmax
        mat[1:(2*n+1),1:(2*n+1),n] = rotmat(Btns,n,nfs[n])
    end
    return mat
end

end
